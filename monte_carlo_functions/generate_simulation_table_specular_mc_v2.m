function generate_simulation_table_specular_mc_v2(sim_table_params, scenarios_file)
    % GENERATE_SIMULATION_TABLE_SPECULAR_MC_V2
    % Generates Monte Carlo scenarios using pure SPECULAR IMAGE SOURCE physics.
    %
    % UPDATED: 
    % - Source Angles are now snapped to the 'omega_grid' (BSM Grid).
    % - Dynamically handles "Inside" and "Outside" FoV conditions.
    % - Dynamically saves to the designated scenarios_file.
    
    % --- CONFIGURATION ---
    N_mc_per_config = sim_table_params.N_mc_per_config;
    base_path = sim_table_params.base_path;
    timit_root = sim_table_params.timit_root;
    diffuese_flag = sim_table_params.diffuese_flag;
    specular_flag = sim_table_params.specular_flag;
    fov_condition = sim_table_params.fov_condition;

    drr_range = sim_table_params.drr_range;
    snr_range = sim_table_params.snr_range;
    
    % Flag for Infinity DRR (Direct Only) to bypass strict DRR matching
    is_inf_drr = (drr_range(1) > 50); 

    min_speech_duration = 3.5; 
    target_counts = sim_table_params.target_counts;
    room_types = sim_table_params.room_types;
     
    % Constraints
    min_dist_from_array = 0.25; 
    wall_margin = 0.5; 
    min_inter_source_dist = 0.25;
    
    % Grid Setup
    source_distribution = 4;
    Q = 36;
    sector_ang = [];
    [el_grid, az_grid] = BSM_toolbox.BSMgrid(source_distribution, Q, sector_ang);
    omega_grid = [el_grid, az_grid]; 
    
    % --- CALIBRATION FACTORS ---
    drr_offset_db = 1.0;  
    snr_offset_db = 2.3;  
    
    % Build Library
    timit_library = get_valid_timit_files(timit_root, min_speech_duration);
    if isempty(timit_library), error('No valid TIMIT files found!'); end
    scenarios = [];
    cnt = 1;
    
    fprintf('Generating Scenarios (FoV: %s, Grid-Based Angles)...\n', fov_condition);
    
    for r_idx = 1:length(room_types)
        curr_room_type = room_types(r_idx);
        
        % Get Room Physics
        roomsetup = setup_room_params(char(curr_room_type),diffuese_flag,specular_flag);
        roomDim   = roomsetup.roomDim;
        arrayPos  = roomsetup.arrayPos;
        
        % Alpha from Reflection Amplitude
        alpha_wall = 1 - (roomsetup.R).^2; 
        
        % Babble Position
        babble_pos = [roomDim(1)-wall_margin, roomDim(2)-wall_margin, 1.7]; 
        
        for t_idx = 1:length(target_counts)
            n_targets = target_counts(t_idx);
            
            for mc = 1:N_mc_per_config
                
                % 1. GENERATE POSITIONS (Grid Selection + DRR Rejection)
                target_positions = zeros(n_targets, 3);
                valid_scene = false; 
                scene_attempt = 0;
                
                % Attempt to build a full valid scene
                while ~valid_scene && scene_attempt < 100
                    targets_found = 0;
                    
                    % Find targets one by one
                    for t = 1:n_targets
                        valid_target = false;
                        target_attempt = 0;
                        
                        % Try to place the current single target
                        while ~valid_target && target_attempt < 1000
                            cand = generate_candidates_fov(1, roomDim, arrayPos, ...
                                                           min_dist_from_array, wall_margin, omega_grid, fov_condition);
                            
                            % Only check distance against targets we have ALREADY found
                            existing = target_positions(1:targets_found, :);
                            
                            if check_single_constraint(cand, existing, babble_pos, arrayPos, min_dist_from_array, min_inter_source_dist)
                                % Check DRR Constraint
                                val_drr_theory = calculate_specular_drr(roomDim, arrayPos, cand, alpha_wall);
                                val_drr_predicted = val_drr_theory + drr_offset_db;
                                
                                % Accept if it's within range OR if we are doing Infinity DRR
                                if is_inf_drr || (val_drr_predicted >= drr_range(1) && val_drr_predicted <= drr_range(2))
                                    valid_target = true;
                                    target_positions(t, :) = cand;
                                    targets_found = targets_found + 1;
                                end
                            end
                            target_attempt = target_attempt + 1;
                        end
                        
                        if ~valid_target
                            break; 
                        end
                    end
                    
                    if targets_found == n_targets
                        valid_scene = true;
                    end
                    scene_attempt = scene_attempt + 1;
                end
                
                if ~valid_scene
                    fprintf('Failed to find valid positions for %d targets after maximum attempts. Skipping MC.\n', n_targets);
                    continue; 
                end
                
                % 2. AUDIO GAIN SETUP
                target_gains = 0.8 + (1.2-0.8).*rand(n_targets, 1); 
                
                % 3. CALCULATE POWER (SPECULAR)
                total_tgt_power = 0;
                for t = 1:n_targets
                    e_tgt = get_specular_energy(roomDim, arrayPos, target_positions(t,:), alpha_wall);
                    total_tgt_power = total_tgt_power + (target_gains(t)^2 * e_tgt);
                end
                
                e_babble = get_specular_energy(roomDim, arrayPos, babble_pos, alpha_wall);
                
                % 4. SOLVE FOR BABBLE GAIN (WITH SNR CALIBRATION)
                desired_snr = snr_range(1) + (snr_range(2)-snr_range(1)) * rand();
                
                req_babble_power = total_tgt_power / 10^(desired_snr/10);
                req_babble_power = req_babble_power * 10^(snr_offset_db/10); % Apply Offset
                
                babble_gain = sqrt( req_babble_power / e_babble );
                
                % Save
                selected_files = timit_library(randsample(length(timit_library), n_targets));
                
                % Ensure a clean output value for DRR if it was bypassed for Infinity
                if is_inf_drr, val_drr_predicted = Inf; end
                
                newRow = struct('ID', cnt, 'RoomType', curr_room_type, 'NumTargets', n_targets, ...
                    'BabblePos', babble_pos, 'BabbleGain', babble_gain, ...
                    'TargetPositions', target_positions, 'TargetGains', target_gains, ...
                    'AudioFiles', string(selected_files), ...
                    'DesiredSNR', desired_snr, ...
                    'DesiredDRR', val_drr_predicted); 
                
                scenarios = [scenarios; newRow]; %#ok<AGROW>
                cnt = cnt + 1;
            end
        end
    end
    
    scenarioTable = struct2table(scenarios);
    save(scenarios_file, 'scenarioTable');
    fprintf('Generated %d scenarios saved to: %s\n', height(scenarioTable), scenarios_file);
end

% --- PHYSICS HELPER: SPECULAR DRR ---
function drr_val = calculate_specular_drr(L, mic, src, alpha)
    d_sq = sum((mic - src).^2);
    e_direct = 1/d_sq;
    e_total = get_specular_energy(L, mic, src, alpha);
    e_reverb = e_total - e_direct;
    if e_reverb < 1e-15, e_reverb = 1e-15; end
    drr_val = 10 * log10(e_direct / e_reverb);
end

% --- PHYSICS HELPER: SPECULAR ENERGY ---
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

% --- POSITION GENERATOR (FOV DYNAMIC) ---
function candidates = generate_candidates_fov(n, roomDim, arrayPos, min_dist, margin, omega_grid, fov_cond)
    % Generates 'n' positions locked to omega_grid based on fov_cond ("Inside" or "Outside")
    
    candidates = zeros(n, 3);
    deg2rad = pi/180;
    
    % 1. Define Constraints
    az_limit_val = 45 * deg2rad; 
    el_center = 90 * deg2rad; 
    el_dev    = 10 * deg2rad;
    el_limits = [el_center - el_dev, el_center + el_dev];
    
    el_vals = omega_grid(:, 1);
    az_vals = omega_grid(:, 2);
    
    % Wrap azimuth values to [-pi, pi] to handle potential 0-360 degree outputs
    az_vals_wrapped = wrapToPi(az_vals);
    
    % 2. Filter Grid Points based on FoV Condition
    if strcmpi(fov_cond, 'Inside')
        % Frontal Cone: between -45 and +45
        valid_indices = find(az_vals_wrapped >= -az_limit_val & az_vals_wrapped <= az_limit_val & ...
                             el_vals >= el_limits(1) & el_vals <= el_limits(2));
    else
        % Outside FoV: strictly outside -45 to +45
        valid_indices = find((az_vals_wrapped < -az_limit_val | az_vals_wrapped > az_limit_val) & ...
                             el_vals >= el_limits(1) & el_vals <= el_limits(2));
    end
    
    if isempty(valid_indices)
        error('Grid error: No valid points found for FoV condition "%s".', fov_cond);
    end

    max_room_dist = norm(roomDim);

    for i = 1:n
        % 3. Pick Random Valid Direction from Grid
        idx = valid_indices(randi(length(valid_indices)));
        el = el_vals(idx);
        az = az_vals_wrapped(idx);
        
        % 4. Pick Random Distance
        r = min_dist + (max_room_dist - min_dist) * rand();
        
        % 5. Convert to Cartesian
        dx = r * sin(el) * cos(az);
        dy = r * sin(el) * sin(az);
        dz = r * cos(el);
        
        cand = arrayPos + [dx, dy, dz];
        
        % 6. Clamp to Room
        cand(1) = max(margin, min(roomDim(1)-margin, cand(1)));
        cand(2) = max(margin, min(roomDim(2)-margin, cand(2)));
        cand(3) = max(1.0,    min(1.8,               cand(3))); 
        
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

function ok = check_single_constraint(cand, existing_targets, babble_pos, arrayPos, d_arr, d_src)
    ok = true;
    if norm(cand - arrayPos) < d_arr
        ok = false; return;
    end
    if norm(cand - babble_pos) < d_src
        ok = false; return;
    end
    for i = 1:size(existing_targets, 1)
        if norm(cand - existing_targets(i, :)) < d_src
            ok = false; return;
        end
    end
end