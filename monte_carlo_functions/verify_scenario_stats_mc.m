function verify_scenario_stats_mc(base_path)
    % VERIFY_SCENARIO_STATS_MC
    % Verifies Monte Carlo scenarios generated using specular image source physics.
    
    load(fullfile(base_path, 'simulation_scenarios.mat'), 'scenarioTable');
    stats = struct('Room', [], 'SNR', [], 'DRR', [], 'Dist', []);
    count = 1;
    
    % Calibration Offsets (Must match generator)
    drr_offset_db = 1.0;  % Reduced from 2.5
    snr_offset_db = 2.3;  % Reduced from 4.0
    
    for i = 1:height(scenarioTable)
        row = scenarioTable(i,:);
        
        % Handle table cells dynamically if stored as such
        if iscell(row.RoomType), room_type = char(row.RoomType{1}); else, room_type = char(row.RoomType); end
        if iscell(row.BabblePos), b_pos = row.BabblePos{1}; else, b_pos = row.BabblePos; end
        if iscell(row.TargetPositions), t_mat = row.TargetPositions{1}; else, t_mat = row.TargetPositions; end
        if iscell(row.TargetGains), g_vec = row.TargetGains{1}; else, g_vec = row.TargetGains; end
        
        roomsetup = setup_room_params(room_type,true,true);
        roomDim = roomsetup.roomDim;
        arrayPos = roomsetup.arrayPos;
        
        % Calculate Alpha from Reflection Amplitude
        alpha_wall = 1 - (roomsetup.R).^2; 
        
        % 1. Calculate Actual Babble Power (Physics Aware)
        e_babble = get_specular_energy(roomDim, arrayPos, b_pos, alpha_wall);
        actual_babble_power = (row.BabbleGain^2) * e_babble;
        
        % 2. Calculate Total Target Power (Physics Aware)
        total_tgt_power = 0;
        for t = 1:row.NumTargets
            e_tgt = get_specular_energy(roomDim, arrayPos, t_mat(t,:), alpha_wall);
            total_tgt_power = total_tgt_power + (g_vec(t)^2 * e_tgt);
        end
        
        % 3. Calculate Scene SNR (Theory + Calibration Offset)
        snr_theory = 10 * log10(total_tgt_power / actual_babble_power);
        snr_predicted = snr_theory + snr_offset_db;
        
        % 4. Calculate DRR per Target
        for t = 1:row.NumTargets
            d_target = norm(t_mat(t,:) - arrayPos);
            
            % DRR (Theory + Calibration Offset)
            val_drr_theory = calculate_specular_drr(roomDim, arrayPos, t_mat(t,:), alpha_wall);
            val_drr_predicted = val_drr_theory + drr_offset_db;
            
            % Save to Stats
            stats(count).Room = string(room_type);
            stats(count).SNR  = snr_predicted;
            stats(count).DRR  = val_drr_predicted;
            stats(count).Dist = d_target;
            count = count + 1;
        end
    end
    
    T = struct2table(stats);
    
    % --- PLOTTING ---
    figure('Name', 'Scenario Statistics (Specular MC)');
    room_order = ["small", "medium", "large"];
    cols = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125];
    
    subplot(1,3,1); hold on;
    for i=1:3
        idx = T.Room == room_order(i);
        if any(idx)
            scatter(T.Dist(idx), T.DRR(idx), 15, cols(i,:), 'filled', 'DisplayName', room_order(i));
        end
    end
    yline([-10, 10], '--k', 'Alpha', 0.3, 'HandleVisibility', 'off');
    xlabel('Distance (m)'); ylabel('DRR (dB)'); title('Target DRR vs Distance'); grid on; legend;
    
    subplot(1,3,2); hold on;
    for i=1:3
        idx = T.Room == room_order(i);
        if any(idx)
            histogram(T.SNR(idx), 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', cols(i,:), 'DisplayName', room_order(i));
        end
    end
    xlabel('SNR (dB)'); title('Global Scene SNR'); grid on; legend;

    subplot(1,3,3); hold on;
    for i=1:3
        idx = T.Room == room_order(i);
        if any(idx)
            histogram(T.DRR(idx), 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', cols(i,:), 'DisplayName', room_order(i));
        end
    end
    xlabel('DRR (dB)'); title('Target DRR Distribution'); grid on; legend;
end

% =========================================================================
% PHYSICS HELPERS (MIRRORED FROM GENERATOR SCRIPT)
% =========================================================================

function drr_val = calculate_specular_drr(L, mic, src, alpha)
    d_sq = sum((mic - src).^2);
    e_direct = 1/d_sq;
    e_total = get_specular_energy(L, mic, src, alpha);
    e_reverb = e_total - e_direct;
    if e_reverb < 1e-15, e_reverb = 1e-15; end
    drr_val = 10 * log10(e_direct / e_reverb);
end

function E_total = get_specular_energy(L, mic, src, alpha)
    order = 3; 
    R_energy = 1 - alpha;
    d_sq = sum((mic - src).^2);
    E_total = 1/d_sq; 
    for x = -order:order
        for y = -order:order
            for z = -order:order
                if x==0 && y==0 && z==0, continue; end 
                sx = get_img_coord(src(1), L(1), x);
                sy = get_img_coord(src(2), L(2), y);
                sz = get_img_coord(src(3), L(3), z);
                dist_sq = (mic(1)-sx)^2 + (mic(2)-sy)^2 + (mic(3)-sz)^2;
                N_refl = abs(x) + abs(y) + abs(z);
                E_total = E_total + (R_energy^N_refl) / dist_sq;
            end
        end
    end
end

function coord = get_img_coord(s, L, n)
    if mod(n, 2) == 0, coord = n*L + s; else, coord = n*L + (L - s); end
end