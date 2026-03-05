function generate_simulation_table_simp(n_distance_steps, base_path, timit_root)
    % MODIFIED: Generates a deterministic "Distance Sweep" along the frontal axis.
    % INCLUDES CALIBRATION for "Tom" Simulator (Image Method w/o Scattering)
    
    % Configuration
    min_speech_duration = 3.5; 
    room_types = ["small", "medium", "large"];
    
    % Forced simplification for distance sweep analysis
    target_counts = [1];
    
    % Constraints
    min_dist_from_array = 0.3; 
    wall_margin = 0.5; 
    
    % SNR Settings
    snr_fixed = 17.5; % dB
    
    % --- CALIBRATION PARAMETERS (Derived from Validation Plots) ---
    % Problem: "Tom" sim produces less reverb than diffuse theory.
    % Fix 1: Effective rc is larger (DRR is ~3.5dB higher -> factor 1.5)
    rc_calibration_factor = 1.5; 
    
    % Fix 2: Even with rc fix, Babble (corner) is too quiet. 
    % SNR Error was ~7dB. The rc_factor fixes ~3.5dB of it. 
    % We need an extra boost to force SNR alignment.
    babble_gain_tweak_db  = 3.5; 
    
    % Build Library
    timit_library = get_valid_timit_files(timit_root, min_speech_duration);
    if isempty(timit_library), error('No valid TIMIT files found!'); end
    scenarios = [];
    cnt = 1;
    
    fprintf('Generating Distance Sweep Scenarios (Calibrated for Image Method)...\n');
    
    for r_idx = 1:length(room_types)
        curr_room_type = room_types(r_idx);
        
        % Get Room Physics
        roomsetup = setup_room_params(char(curr_room_type));
        roomDim  = roomsetup.roomDim;
        arrayPos = roomsetup.arrayPos;
        
        % Apply Calibration to Critical Distance
        % This acknowledges that the room "sounds" drier than Sabine theory
        r_c_theoretical = roomsetup.critical_distance;
        r_c_effective   = r_c_theoretical * rc_calibration_factor;
        
        % Fixed Babble Position (Corner)
        babble_pos = [roomDim(1)-wall_margin, roomDim(2)-wall_margin, 1.7]; 
        d_babble   = norm(babble_pos - arrayPos);
        
        % Distance Steps Setup
        dir_vec = [1, 0, 0]; 
        dist_to_wall = roomDim(1) - arrayPos(1);
        max_dist = dist_to_wall - wall_margin;
        
        % Sweep from min_dist to effective critical distance (or max_dist)
        dist_steps = linspace(min_dist_from_array, min(max_dist, r_c_effective*1.5), n_distance_steps);
        
        for t_idx = 1:length(target_counts)
            n_targets = target_counts(t_idx);
            
            for k = 1:n_distance_steps
                curr_dist = dist_steps(k);
                
                % 1. Set Position
                target_positions = arrayPos + curr_dist * dir_vec;
                
                % 2. Audio Gain
                target_gains = 1.0; 
                
                % 3. Physics & SNR (Using Effective Parameters)
                % Power ~ Gain^2 * (1/d^2 + 1/rc_eff^2)
                tgt_power = target_gains^2 * ( (1/curr_dist^2) + (1/r_c_effective^2) );
                
                % Solve for Required Babble Power
                req_babble_power = tgt_power / 10^(snr_fixed/10);
                
                % Apply Extra Heuristic Boost to Babble Power (to fix the +7dB SNR offset)
                req_babble_power = req_babble_power * 10^(babble_gain_tweak_db/10);
                
                % Calculate Babble Gain
                babble_geometric_factor = (1/d_babble^2) + (1/r_c_effective^2);
                babble_gain = sqrt( req_babble_power / babble_geometric_factor );
                
                selected_files = timit_library(randsample(length(timit_library), n_targets));
                
                newRow = struct('ID', cnt, 'RoomType', curr_room_type, 'NumTargets', n_targets, ...
                    'BabblePos', babble_pos, 'BabbleGain', babble_gain, ...
                    'TargetPositions', target_positions, 'TargetGains', target_gains, ...
                    'AudioFiles', string(selected_files));
                
                scenarios = [scenarios; newRow]; %#ok<AGROW>
                cnt = cnt + 1;
            end
        end
    end
    
    scenarioTable = struct2table(scenarios);
    save(fullfile(base_path, 'simulation_scenarios.mat'), 'scenarioTable');
    fprintf('Generated %d scenarios.\n', height(scenarioTable));
end

% --- HELPERS ---
function file_list = get_valid_timit_files(root_path, min_dur)
    cache_file = fullfile(root_path, 'valid_file_cache.mat');
    if exist(cache_file, 'file')
        load(cache_file, 'file_list', 'cached_dur');
        if cached_dur == min_dur, return; end
    end
    files = dir(fullfile(root_path, '**/*.wav'));
    valid = {}; c=0;
    for i=1:length(files)
        try, if audioinfo(fullfile(files(i).folder, files(i).name)).Duration >= min_dur, c=c+1; valid{c,1}=fullfile(files(i).folder, files(i).name); end, catch, end %#ok<AGROW>
    end
    file_list = string(valid); cached_dur = min_dur;
    save(cache_file, 'file_list', 'cached_dur');
end