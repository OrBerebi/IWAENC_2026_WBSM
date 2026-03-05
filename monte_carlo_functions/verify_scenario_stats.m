function verify_scenario_stats(base_path)
    load(fullfile(base_path, 'simulation_scenarios.mat'), 'scenarioTable');
    stats = struct('Room', [], 'SNR', [], 'DRR', [], 'Dist', []);
    count = 1;
    
    for i = 1:height(scenarioTable)
        row = scenarioTable(i,:);
        roomsetup = setup_room_params(char(row.RoomType));
        r_c = roomsetup.critical_distance; arrayPos = roomsetup.arrayPos;
        
        % Babble Power (Physics Aware)
        if iscell(row.BabblePos), b_pos = row.BabblePos{1}; else, b_pos = row.BabblePos; end
        d_babble = norm(b_pos - arrayPos);
        % Power ~ Gain^2 * (1/d^2 + 1/rc^2)
        pow_babble = (row.BabbleGain^2) * ( (1/d_babble^2) + (1/r_c^2) );
        
        if iscell(row.TargetPositions), t_mat = row.TargetPositions{1}; else, t_mat = row.TargetPositions; end
        if iscell(row.TargetGains), g_vec = row.TargetGains{1}; else, g_vec = row.TargetGains; end
        
        for t = 1:row.NumTargets
            d_target = norm(t_mat(t,:) - arrayPos);
            
            % Target Power (Physics Aware)
            pow_target = (g_vec(t)^2) * ( (1/d_target^2) + (1/r_c^2) );
            
            % Metrics
            snr_val = 10 * log10(pow_target / pow_babble); % Ratio of Powers = 10log10
            drr_val = 20 * log10(r_c / d_target);
            
            stats(count).Room = string(row.RoomType);
            stats(count).SNR  = snr_val;
            stats(count).DRR  = drr_val;
            stats(count).Dist = d_target;
            count = count + 1;
        end
    end
    T = struct2table(stats);
    
    % Plotting
    figure('Name', 'Scenario Statistics');
    room_order = ["small", "medium", "large"];
    cols = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125];
    
    subplot(1,3,1); hold on;
    for i=1:3, idx=T.Room==room_order(i); scatter(T.Dist(idx), T.DRR(idx), 15, cols(i,:), 'filled', 'DisplayName', room_order(i)); end
    yline([-10, 10], '--k', 'Alpha', 0.3, 'HandleVisibility', 'off');
    xlabel('Distance (m)'); ylabel('DRR (dB)'); title('DRR vs Distance'); grid on; legend;
    
    subplot(1,3,2); hold on;
    for i=1:3, histogram(T.SNR(T.Room==room_order(i)), 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', cols(i,:), 'DisplayName', room_order(i)); end
    xlabel('SNR (dB)'); title('SNR'); grid on; legend;

    subplot(1,3,3); hold on;
    for i=1:3, histogram(T.DRR(T.Room==room_order(i)), 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', cols(i,:), 'DisplayName', room_order(i)); end
    xlabel('DRR (dB)'); title('DRR'); grid on; legend;
end