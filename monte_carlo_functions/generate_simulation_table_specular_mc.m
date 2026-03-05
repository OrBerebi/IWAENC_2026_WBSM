function generate_simulation_table_specular_mc(N_mc_per_config, base_path, timit_root)
    % GENERATE_SIMULATION_TABLE_SPECULAR_MC
    % Generates Monte Carlo scenarios using pure SPECULAR IMAGE SOURCE physics.
    %
    % 1. Positioning: Frontal Cone (+/-45 deg).
    % 2. Distance: Determined by calculating the EXACT SPECULAR DRR (Order 3)
    %    and accepting only positions within 'drr_range'.
    % 3. Gains: Calculated using Specular Energy summation.
    
    % --- CONFIGURATION ---
    min_speech_duration = 3.5; 
    room_types = ["small", "medium", "large"];
    target_counts = [1, 2, 3]; 
    
    % Constraints
    min_dist_from_array = 0.3; 
    wall_margin = 0.5; 
    min_inter_source_dist = 0.25;
    
    % RANGES
    drr_range = [8, 20];    % dB (Selection Criteria based on Specular Energy)
    snr_range = [15, 20];   % dB (Target-to-Babble Ratio)
    
    % Build Library
    timit_library = get_valid_timit_files(timit_root, min_speech_duration);
    if isempty(timit_library), error('No valid TIMIT files found!'); end
    scenarios = [];
    cnt = 1;
    
    fprintf('Generating Scenarios (Specular DRR Selection, Frontal +/-45)...\n');
    
    for r_idx = 1:length(room_types)
        curr_room_type = room_types(r_idx);
        
        % Get Room Physics
        roomsetup = setup_room_params(char(curr_room_type));
        roomDim   = roomsetup.roomDim;
        arrayPos  = roomsetup.arrayPos;
        
        % EXTRACT ALPHA (ABSORPTION) FROM AMPLITUDE
        % roomsetup.R is Reflection Amplitude. alpha = 1 - |r|^2
        alpha_wall = 1 - (roomsetup.R).^2; 
        
        % Babble Position (Fixed in Corner)
        babble_pos = [roomDim(1)-wall_margin, roomDim(2)-wall_margin, 1.7]; 
        
        for t_idx = 1:length(target_counts)
            n_targets = target_counts(t_idx);
            
            for mc = 1:N_mc_per_config
                
                % 1. GENERATE POSITIONS (Rejection Sampling based on Specular DRR)
                target_positions = zeros(n_targets, 3);
                valid_scene = false; attempt = 0;
                
                while ~valid_scene && attempt < 5000
                    % Generate random candidates in Frontal Cone
                    candidates = generate_frontal_candidates(n_targets, roomDim, arrayPos, ...
                                                             min_dist_from_array, wall_margin);
                    
                    % Check Physical Constraints (Walls, Inter-source)
                    if check_constraints(candidates, babble_pos, arrayPos, min_dist_from_array, min_inter_source_dist)
                        
                        % Check DRR Constraint (The Key Step)
                        drr_ok = true;
                        for t = 1:n_targets
                            % Calculate Exact Specular DRR for this candidate
                            val_drr = calculate_specular_drr(roomDim, arrayPos, candidates(t,:), alpha_wall);
                            
                            if val_drr < drr_range(1) || val_drr > drr_range(2)
                                drr_ok = false; break;
                            end
                        end
                        
                        if drr_ok
                            target_positions = candidates;
                            valid_scene = true;
                        end
                    end
                    attempt = attempt + 1;
                end
                
                if ~valid_scene
                    % fprintf('Warning: Could not find valid positions for Scene %d (Room %s)\n', cnt, curr_room_type);
                    continue; 
                end
                
                % 2. AUDIO GAIN SETUP
                target_gains = 0.8 + (1.2-0.8).*rand(n_targets, 1); 
                
                % 3. CALCULATE TARGET POWER (SPECULAR)
                total_tgt_power = 0;
                for t = 1:n_targets
                    e_tgt = get_specular_energy(roomDim, arrayPos, target_positions(t,:), alpha_wall);
                    total_tgt_power = total_tgt_power + (target_gains(t)^2 * e_tgt);
                end
                
                % 4. CALCULATE BABBLE POWER (SPECULAR)
                e_babble = get_specular_energy(roomDim, arrayPos, babble_pos, alpha_wall);
                
                % 5. SOLVE FOR BABBLE GAIN
                desired_snr = snr_range(1) + (snr_range(2)-snr_range(1)) * rand();
                
                % P_babble_req = P_tgt_total / 10^(SNR/10)
                req_babble_power = total_tgt_power / 10^(desired_snr/10);
                babble_gain = sqrt( req_babble_power / e_babble );
                
                % Save
                selected_files = timit_library(randsample(length(timit_library), n_targets));
                
                newRow = struct('ID', cnt, 'RoomType', curr_room_type, 'NumTargets', n_targets, ...
                    'BabblePos', babble_pos, 'BabbleGain', babble_gain, ...
                    'TargetPositions', target_positions, 'TargetGains', target_gains, ...
                    'AudioFiles', string(selected_files), ...
                    'DesiredSNR', desired_snr); 
                
                scenarios = [scenarios; newRow]; %#ok<AGROW>
                cnt = cnt + 1;
            end
        end
    end
    
    scenarioTable = struct2table(scenarios);
    save(fullfile(base_path, 'simulation_scenarios.mat'), 'scenarioTable');
    fprintf('Generated %d scenarios.\n', height(scenarioTable));
end

% --- PHYSICS HELPER: SPECULAR DRR ---
function drr_val = calculate_specular_drr(L, mic, src, alpha)
    % Calculates DRR using Deterministic Image Sources (Order 3)
    % DRR = 10 * log10( E_direct / E_reverb )
    
    % Direct Energy
    d_sq = sum((mic - src).^2);
    e_direct = 1/d_sq;
    
    % Total Specular Energy
    e_total = get_specular_energy(L, mic, src, alpha);
    
    % Reverb Energy
    e_reverb = e_total - e_direct;
    if e_reverb < 1e-15, e_reverb = 1e-15; end
    
    drr_val = 10 * log10(e_direct / e_reverb);
end

% --- PHYSICS HELPER: SPECULAR ENERGY ---
function E_total = get_specular_energy(L, mic, src, alpha)
    % Calculates Total Energy (Direct + Order 3 Reflections)
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

% --- POSITION GENERATOR ---
function candidates = generate_frontal_candidates(n, roomDim, arrayPos, min_dist, margin)
    % Generates 'n' random positions in Frontal Cone
    % Does NOT use DRR to set distance (that is checked later).
    % Distance is sampled uniformly to cover the room.
    
    candidates = zeros(n, 3);
    deg2rad = pi/180;
    
    % Constraints
    az_limits = [-45, 45] * deg2rad; 
    el_center = 90 * deg2rad; 
    el_dev    = 10 * deg2rad;
    el_limits = [el_center - el_dev, el_center + el_dev];
    
    % Distance Sampling Max (Corner of room)
    max_room_dist = norm(roomDim);

    for i = 1:n
        % Random Distance (Uniform sampling usually sufficient for MC)
        r = min_dist + (max_room_dist - min_dist) * rand();
        
        az = az_limits(1) + (az_limits(2) - az_limits(1)) * rand();
        el = el_limits(1) + (el_limits(2) - el_limits(1)) * rand();
        
        dx = r * sin(el) * cos(az);
        dy = r * sin(el) * sin(az);
        dz = r * cos(el);
        
        cand = arrayPos + [dx, dy, dz];
        
        % Clamp to room dimensions (if outside, logic in main loop will regenerate or just cap)
        % Simple clamp here for candidate validity:
        cand(1) = max(margin, min(roomDim(1)-margin, cand(1)));
        cand(2) = max(margin, min(roomDim(2)-margin, cand(2)));
        cand(3) = max(1.0,    min(1.8,               cand(3))); % Height constraint
        
        candidates(i,:) = cand;
    end
end

% --- UTILS ---
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