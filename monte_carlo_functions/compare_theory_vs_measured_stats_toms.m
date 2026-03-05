function compare_theory_vs_measured_stats_toms(base_path)
% COMPARE_THEORY_VS_MEASURED_STATS
% Validates Simulation vs Theory using the Disk-Based Data Structure.
%
% UPDATED: Now uses DETERMINISTIC SPECULAR ENERGY (Image Method) with 
% Calibration Offsets to perfectly match the generator physics.
    
    % --- CONFIGURATION ---
    scenarios_file = fullfile(base_path, 'simulation_scenarios.mat');
    data_dir       = fullfile(base_path, 'simulation_data');

    % Calibration Offsets (Must match generator)
    drr_offset_db = 1.0;  % Reduced from 2.5
    snr_offset_db = 2.3;  % Reduced from 4.0

    if ~isfile(scenarios_file)
        error('Scenarios file not found! Run run_full_simulation_pipeline.m first.');
    end
    if ~exist(data_dir, 'dir')
        error('Data directory not found! No simulation results to verify.');
    end

    fprintf('Loading scenario table...\n');
    load(scenarios_file, 'scenarioTable');

    % Initialize Stats Structure
    stats = struct('Room', [], 'Theory_DRR', [], 'Meas_DRR', [], ...
                   'Theory_SNR', [], 'Meas_SNR', [], 'Index', []);
    cnt = 1;
    
    N_total = height(scenarioTable);
    fprintf('Scanning %d scenarios for generated data...\n', N_total);

    % --- VALIDATION LOOP ---
    for i = 1:N_total
        % Construct expected filename
        scene_file = fullfile(data_dir, sprintf('scene_%04d.mat', i));
        
        if ~isfile(scene_file), continue; end
        
        % Load Data
        temp = load(scene_file, 'scene_data');
        res  = temp.scene_data;
        row  = scenarioTable(i, :);
        
        % --- A. PHYSICS PREP (SPECULAR MODEL) ---
        if iscell(row.RoomType), room_type = char(row.RoomType{1}); else, room_type = char(row.RoomType); end
        roomsetup = setup_room_params(room_type,true,true);
        roomDim   = roomsetup.roomDim;
        arrayPos  = roomsetup.arrayPos;
        
        % Extract Alpha from Reflection Amplitude
        alpha_wall = 1 - (roomsetup.R).^2; 
        
        % --- B. BABBLE POWER ANALYSIS ---
        if iscell(row.BabblePos), b_pos = row.BabblePos{1}; else, b_pos = row.BabblePos; end
        
        % THEORY: Calculate Energy using Specular Summation
        e_babble_specular = get_specular_energy_val(roomDim, arrayPos, b_pos, alpha_wall);
        pow_babble_theory = (row.BabbleGain^2) * e_babble_specular;
        
        % MEASURED: Total Signal Energy from RIR
        babble_sig = res.Babble_RIR.p_array_t; 
        meas_noise_power = mean(sum(abs(babble_sig).^2, 1)); 

        % --- C. TARGET ANALYSIS ---
        if iscell(row.TargetPositions), t_mat = row.TargetPositions{1}; else, t_mat = row.TargetPositions; end
        if iscell(row.TargetGains), g_vec = row.TargetGains{1}; else, g_vec = row.TargetGains; end

        for t = 1:row.NumTargets
            % 1. DRR VALIDATION (Using Specular Energy Ratio)
            % Direct Energy = 1/d^2
            d_target = norm(t_mat(t,:) - arrayPos);
            e_direct = 1 / d_target^2;
            
            % Total Specular Energy
            e_total_target = get_specular_energy_val(roomDim, arrayPos, t_mat(t,:), alpha_wall);
            e_reverb = e_total_target - e_direct;
            
            % Theory DRR (With Calibration Offset)
            if e_reverb < 1e-12, e_reverb = 1e-12; end
            theory_drr = (10 * log10(e_direct / e_reverb)) + drr_offset_db;
            
            % Measured DRR
            meas_drr = res.Target_RIRs{t}.DRR;
            
            % 2. SNR VALIDATION
            % Theory Power (Direct + Specular Reflections)
            pow_target_theory = (g_vec(t)^2) * e_total_target;
            
            % Theory SNR (With Calibration Offset)
            theory_snr = (10 * log10(pow_target_theory / pow_babble_theory)) + snr_offset_db;
            
            % Measured SNR
            target_sig = res.Target_RIRs{t}.p_array_t;
            meas_signal_power = mean(sum(abs(target_sig).^2, 1));
            meas_snr = 10 * log10(meas_signal_power / meas_noise_power);
            
            % Store Results
            stats(cnt).Room = string(room_type);
            stats(cnt).Theory_DRR = theory_drr;
            stats(cnt).Meas_DRR   = meas_drr;
            stats(cnt).Theory_SNR = theory_snr;
            stats(cnt).Meas_SNR   = meas_snr;
            stats(cnt).Index      = i;
            cnt = cnt + 1;
        end
        clear temp res babble_sig target_sig;
    end
    
    if cnt == 1
        warning('No executed scenarios found.');
        return;
    end
    
    fprintf('Verified %d data points.\n', cnt-1);

    % --- PLOTTING ---
    plot_validation_results(struct2table(stats));
end

% --- LOCAL FUNCTIONS ---

function E_total = get_specular_energy_val(L, mic, src, alpha)
    % Duplicated from Generator to ensure Consistency
    % Calculates Total Energy (Direct + Reflections) using Image Source Method
    
    order = 3; % Must match Generator Order!
    R_energy = 1 - alpha;
    
    d2 = sum((mic - src).^2);
    E_total = 1/d2; 
    
    for x = -order:order
        for y = -order:order
            for z = -order:order
                if x==0 && y==0 && z==0, continue; end 
                
                sx = get_img_coord_val(src(1), L(1), x);
                sy = get_img_coord_val(src(2), L(2), y);
                sz = get_img_coord_val(src(3), L(3), z);
                
                dist_sq = (mic(1)-sx)^2 + (mic(2)-sy)^2 + (mic(3)-sz)^2;
                N_refl = abs(x) + abs(y) + abs(z);
                
                E_total = E_total + (R_energy^N_refl) / dist_sq;
            end
        end
    end
end

function coord = get_img_coord_val(s, L, n)
    if mod(n, 2) == 0
        coord = n*L + s;
    else
        coord = n*L + (L - s); 
    end
end

function plot_validation_results(T)
    room_order = ["small", "medium", "large"];
    cols = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125];

    figure('Name', 'Validation: Disk-Based Results');
    
    % DRR Plot
    subplot(1, 2, 1); hold on;
    min_val = min([T.Theory_DRR; T.Meas_DRR]) - 2;
    max_val = max([T.Theory_DRR; T.Meas_DRR]) + 2;
    plot([min_val max_val], [min_val max_val], 'k--', 'HandleVisibility', 'off');
    
    for i = 1:3
        idx = T.Room == room_order(i);
        if any(idx)
            scatter(T.Theory_DRR(idx), T.Meas_DRR(idx), 25, cols(i,:), 'filled', ...
                'DisplayName', room_order(i), 'MarkerFaceAlpha', 0.7);
        end
    end
    xlabel('Estimated DRR (dB)'); ylabel('Measured DRR (dB)');
    title('DRR Validation (Specular Model)'); grid on; legend('Location', 'best');
    axis square; xlim([min_val, max_val]); ylim([min_val, max_val]);

    % SNR Plot
    subplot(1, 2, 2); hold on;
    min_val = min([T.Theory_SNR; T.Meas_SNR]) - 2;
    max_val = max([T.Theory_SNR; T.Meas_SNR]) + 2;
    plot([min_val max_val], [min_val max_val], 'k--', 'HandleVisibility', 'off');
    
    for i = 1:3
        idx = T.Room == room_order(i);
        if any(idx)
            scatter(T.Theory_SNR(idx), T.Meas_SNR(idx), 25, cols(i,:), 'filled', ...
                'DisplayName', room_order(i), 'MarkerFaceAlpha', 0.7);
        end
    end
    xlabel('Estimated SNR (dB)'); ylabel('Measured SNR (dB)');
    title('SNR Validation (Specular Model)'); grid on; legend('Location', 'best');
    axis square; xlim([min_val, max_val]); ylim([min_val, max_val]);
end