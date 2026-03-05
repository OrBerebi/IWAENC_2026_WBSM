function generate_simulation_table(N_mc_per_config, base_path, timit_root)
    % Configuration
    min_speech_duration = 3.5; 
    room_types = ["small", "medium", "large"];
    target_counts = [1, 2, 3]; 
    
    % Constraints
    min_dist_from_array = 0.3; 
    wall_margin = 0.5; 
    min_inter_source_dist = 0.25;
    drr_range = [10, 15];   % dB
    snr_range = [10, 25];   % dB (Total Energy SNR)
    
    % Build Library
    timit_library = get_valid_timit_files(timit_root, min_speech_duration);
    if isempty(timit_library), error('No valid TIMIT files found!'); end
    scenarios = [];
    cnt = 1;
    
    fprintf('Generating scenarios with FRONTAL constraints (+/-50 az, +/-10 el)...\n');

    for r_idx = 1:length(room_types)
        curr_room_type = room_types(r_idx);
        
        % Get Room Physics
        roomsetup = setup_room_params(char(curr_room_type));
        roomDim  = roomsetup.roomDim;
        arrayPos = roomsetup.arrayPos;
        r_c      = roomsetup.critical_distance;
        
        % Fixed Babble Position
        babble_pos = [roomDim(1)-wall_margin, roomDim(2)-wall_margin, 1.7]; 
        d_babble   = norm(babble_pos - arrayPos);
        
        % DRR Sampling Radii
        r_inner = r_c * 10^(-drr_range(2)/20); % Max DRR
        r_outer = r_c * 10^(-drr_range(1)/20); % Min DRR
        r_inner = max(r_inner, min_dist_from_array);
        r_outer = min(r_outer, norm(roomDim));
        
        for t_idx = 1:length(target_counts)
            n_targets = target_counts(t_idx);
            
            for mc = 1:N_mc_per_config
                % Geometry
                target_positions = zeros(n_targets, 3);
                valid_scene = false; attempt = 0;
                
                while ~valid_scene && attempt < 2000
                    % --- MODIFIED CALL: Generates only in Frontal Cone ---
                    candidates = generate_drr_focused_pos_frontal(n_targets, roomDim, arrayPos, r_inner, r_outer, wall_margin);
                    
                    if check_constraints(candidates, babble_pos, arrayPos, min_dist_from_array, min_inter_source_dist)
                        target_positions = candidates;
                        valid_scene = true;
                    end
                    attempt = attempt + 1;
                end
                
                if ~valid_scene
                    fprintf('Warning: Could not find valid frontal positions for Scene %d (Room: %s)\n', cnt, curr_room_type);
                    continue; 
                end
                
                % Audio Gain Setup
                target_gains = 0.8 + (1.2-0.8).*rand(n_targets, 1); 
                
                % --- PHYSICS-AWARE SNR CALCULATION ---
                tgt_powers = zeros(n_targets, 1);
                for t = 1:n_targets
                    d_t = norm(target_positions(t,:) - arrayPos);
                    tgt_powers(t) = target_gains(t)^2 * ( (1/d_t^2) + (1/r_c^2) );
                end
                avg_tgt_power = mean(tgt_powers);
                
                % Pick Desired SNR
                desired_snr = snr_range(1) + (snr_range(2)-snr_range(1)) * rand();
                
                % Solve for Required Babble Gain
                req_babble_power = avg_tgt_power / 10^(desired_snr/10);
                babble_geometric_factor = (1/d_babble^2) + (1/r_c^2);
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

function pos_list = generate_drr_focused_pos_frontal(n, roomDim, arrayPos, r_min, r_max, margin)
    % MODIFIED: Generates positions only within +/- 50 deg Azimuth and +/- 10 deg Elevation
    pos_list = zeros(n, 3);
    deg2rad = pi/180;
    
    % Frontal Constraints (Assuming 0 deg Azimuth is +X axis)
    az_limits = [-50, 50] * deg2rad; 
    
    % Elevation Constraints (90 deg is horizontal plane)
    el_center = 90 * deg2rad; 
    el_dev    = 10 * deg2rad;
    el_limits = [el_center - el_dev, el_center + el_dev];

    for i = 1:n
        valid = false;
        while ~valid
            % Random Radius
            r = r_min + (r_max - r_min) * rand(); 
            
            % Random Azimuth (Restricted)
            az = az_limits(1) + (az_limits(2) - az_limits(1)) * rand();
            
            % Random Elevation (Restricted)
            el = el_limits(1) + (el_limits(2) - el_limits(1)) * rand();
            
            % Spherical to Cartesian (Front = +X)
            % x = r * sin(el) * cos(az)
            % y = r * sin(el) * sin(az)
            % z = r * cos(el)
            
            dx = r * sin(el) * cos(az);
            dy = r * sin(el) * sin(az);
            dz = r * cos(el);
            
            cand = arrayPos + [dx, dy, dz];
            
            % Boundary Checks (Room Dimensions)
            if cand(1)>margin && cand(1)<roomDim(1)-margin && ...
               cand(2)>margin && cand(2)<roomDim(2)-margin && ...
               cand(3)>1.0 && cand(3)<1.8
                
                pos_list(i,:) = cand; 
                valid = true;
            end
        end
    end
end

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

function ok = check_constraints(t, b, a, d_arr, d_src)
    ok=true; all=[t;b];
    for i=1:size(all,1), if norm(all(i,:)-a)<d_arr, ok=false; return; end, end
    for i=1:size(all,1), for j=i+1:size(all,1), if norm(all(i,:)-all(j,:))<d_src, ok=false; return; end, end, end
end