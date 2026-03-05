function compare_theory_vs_measured_stats(base_path)
% COMPARE_THEORY_VS_MEASURED_STATS
% Validates Simulation vs Theory using the Disk-Based Data Structure.
%
% 1. Scans 'simulation_data' folder for executed scenarios.
% 2. Loads each scenario individually to save RAM.
% 3. Compares Measured SNR/DRR against Physics-Aware Theoretical Models.
% 4. Generates Square Plots for validation.
    
    % --- CONFIGURATION ---
    scenarios_file = fullfile(base_path, 'simulation_scenarios.mat');
    data_dir       = fullfile(base_path, 'simulation_data');

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
        % Construct expected filename (e.g., scene_0042.mat)
        scene_file = fullfile(data_dir, sprintf('scene_%04d.mat', i));
        
        % Skip if this scenario was not run (Random 20% selection)
        if ~isfile(scene_file)
            continue; 
        end
        
        % Load Data (One file at a time)
        temp = load(scene_file, 'scene_data');
        res  = temp.scene_data;
        row  = scenarioTable(i, :);
        
        % --- A. PHYSICS PREP ---
        roomsetup = setup_room_params(char(row.RoomType));
        r_c       = roomsetup.critical_distance;
        arrayPos  = roomsetup.arrayPos;
        
        % --- B. BABBLE POWER ANALYSIS ---
        if iscell(row.BabblePos), b_pos = row.BabblePos{1}; else, b_pos = row.BabblePos; end
        d_babble = norm(b_pos - arrayPos);
        
        % Theory: Direct + Reverb Energy (Physics-Aware Model)
        % Power ~ Gain^2 * (1/d^2 + 1/rc^2)
        pow_babble_theory = (row.BabbleGain^2) * ( (1/d_babble^2) + (1/r_c^2) );
        
        % Measured: Total Signal Energy from RIR
        % Power = Mean Square Amplitude (avg over channels)
        babble_sig = res.Babble_RIR.p_array_t; 
        meas_noise_power = mean(sum(abs(babble_sig).^2, 1)); 

        % --- C. TARGET ANALYSIS ---
        if iscell(row.TargetPositions), t_mat = row.TargetPositions{1}; else, t_mat = row.TargetPositions; end
        if iscell(row.TargetGains), g_vec = row.TargetGains{1}; else, g_vec = row.TargetGains; end

        for t = 1:row.NumTargets
            d_target = norm(t_mat(t,:) - arrayPos);
            
            % 1. DRR VALIDATION
            % Theory (Geometric)
            theory_drr = 20 * log10(r_c / d_target);
            % Measured (From RIR Analysis)
            meas_drr   = res.Target_RIRs{t}.DRR;
            
            % 2. SNR VALIDATION
            % Theory (Ratio of Total Theoretical Energies)
            pow_target_theory = (g_vec(t)^2) * ( (1/d_target^2) + (1/r_c^2) );
            theory_snr = 10 * log10(pow_target_theory / pow_babble_theory);
            
            % Measured (Ratio of Actual RIR Powers)
            target_sig = res.Target_RIRs{t}.p_array_t;
            meas_signal_power = mean(sum(abs(target_sig).^2, 1));
            meas_snr = 10 * log10(meas_signal_power / meas_noise_power);
            
            % Store Results
            stats(cnt).Room = string(row.RoomType);
            stats(cnt).Theory_DRR = theory_drr;
            stats(cnt).Meas_DRR   = meas_drr;
            stats(cnt).Theory_SNR = theory_snr;
            stats(cnt).Meas_SNR   = meas_snr;
            stats(cnt).Index      = i;
            cnt = cnt + 1;
        end
        
        % Clear temporary data to keep RAM low
        clear temp res babble_sig target_sig;
    end
    
    if cnt == 1
        warning('No executed scenarios found in %s', data_dir);
        return;
    end
    
    fprintf('Verified %d data points.\n', cnt-1);

    % --- PLOTTING (SQUARE AXES) ---
    T = struct2table(stats);
    room_order = ["small", "medium", "large"];
    cols = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125];

    figure('Name', 'Validation: Disk-Based Results');
    
    % --- PLOT 1: DRR ---
    subplot(1, 2, 1); hold on;
    min_val = min([T.Theory_DRR; T.Meas_DRR]) - 2;
    max_val = max([T.Theory_DRR; T.Meas_DRR]) + 2;
    
    % Perfect Match Line
    plot([min_val max_val], [min_val max_val], 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    
    for i = 1:3
        idx = T.Room == room_order(i);
        if any(idx)
            scatter(T.Theory_DRR(idx), T.Meas_DRR(idx), 25, cols(i,:), 'filled', ...
                'DisplayName', room_order(i), 'MarkerFaceAlpha', 0.7);
        end
    end
    xlabel('Estimated DRR (dB)'); ylabel('Measured DRR (dB)');
    title('DRR Validation'); grid on; legend('Location', 'best');
    axis square; xlim([min_val, max_val]); ylim([min_val, max_val]);

    % --- PLOT 2: SNR ---
    subplot(1, 2, 2); hold on;
    min_val = min([T.Theory_SNR; T.Meas_SNR]) - 2;
    max_val = max([T.Theory_SNR; T.Meas_SNR]) + 2;
    
    plot([min_val max_val], [min_val max_val], 'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    
    for i = 1:3
        idx = T.Room == room_order(i);
        if any(idx)
            scatter(T.Theory_SNR(idx), T.Meas_SNR(idx), 25, cols(i,:), 'filled', ...
                'DisplayName', room_order(i), 'MarkerFaceAlpha', 0.7);
        end
    end
    xlabel('Estimated SNR (dB)'); ylabel('Measured SNR (dB)');
    title('SNR Validation'); grid on; legend('Location', 'best');
    axis square; xlim([min_val, max_val]); ylim([min_val, max_val]);
    
    sgtitle(sprintf('Validation: Physics vs Simulation (%d Scenarios)', cnt-1));
end