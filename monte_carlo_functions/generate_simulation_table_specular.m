function generate_simulation_table_specular(n_distance_steps, base_path, timit_root)
    % MODIFIED: Generates Distance Sweep using DETERMINISTIC SPECULAR ENERGY.
    % Replaces statistical Diffuse Theory (rc) with Image Source Summation.
    % This matches "Tom" (Image Method w/o Scattering) physics exactly.
    
    % Configuration
    min_speech_duration = 3.5; 
    room_types = ["small", "medium", "large"];
    target_counts = [1];
    
    % Constraints
    min_dist_from_array = 0.3; 
    wall_margin = 0.5; 
    
    % Fixed SNR
    snr_fixed = 17.5; % dB
    
    % Build Library
    timit_library = get_valid_timit_files(timit_root, min_speech_duration);
    if isempty(timit_library), error('No valid TIMIT files found!'); end
    scenarios = [];
    cnt = 1;
    
    fprintf('Generating Scenarios using Specular Image Source Model (Order 3)...\n');
    
    for r_idx = 1:length(room_types)
        curr_room_type = room_types(r_idx);
        
        % Get Room Physics
        roomsetup = setup_room_params(char(curr_room_type));
        roomDim   = roomsetup.roomDim;
        arrayPos  = roomsetup.arrayPos;
        
        % EXTRACT ALPHA (ABSORPTION) CORRECTLY
        % roomsetup.R is calculated as (1 - exp(...)) which is Alpha (Absorption).
        alpha_wall = 1- (roomsetup.R).^2; 
        
        % Babble Position
        babble_pos = [roomDim(1)-wall_margin, roomDim(2)-wall_margin, 1.7]; 
        
        % Distance Steps
        dir_vec = [1, 0, 0]; 
        dist_to_wall = roomDim(1) - arrayPos(1);
        max_dist = dist_to_wall - wall_margin;
        
        % Sweep logic (Linear distance steps)
        %dist_steps = linspace(min_dist_from_array, max_dist * 0.8, n_distance_steps);
        dist_steps = linspace(min_dist_from_array, roomsetup.critical_distance, n_distance_steps);
        
        for t_idx = 1:length(target_counts)
            n_targets = target_counts(t_idx);
            
            for k = 1:n_distance_steps
                curr_dist = dist_steps(k);
                
                % 1. Position
                target_positions = arrayPos + curr_dist * dir_vec;
                target_gains = 1.0; 
                
                % 2. CALCULATE TARGET POWER (DIRECT + SPECULAR)
                % We do not use 1/rc^2. We sum the actual reflections.
                tgt_energy = get_specular_energy(roomDim, arrayPos, target_positions, alpha_wall);
                tgt_power_total = target_gains^2 * tgt_energy;
                
                % 3. CALCULATE BABBLE POWER (DIRECT + SPECULAR)
                % Babble is in the corner, so it excites modes differently. 
                % ISM captures this geometry (corner loading) automatically.
                babble_energy = get_specular_energy(roomDim, arrayPos, babble_pos, alpha_wall);
                
                % 4. SOLVE FOR BABBLE GAIN
                % Target_Power / (Babble_Gain^2 * Babble_Energy) = 10^(SNR/10)
                req_babble_power = tgt_power_total / 10^(snr_fixed/10);
                babble_gain = sqrt( req_babble_power / babble_energy );
                
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

% --- SPECULAR ENERGY MODEL ---
function E_total = get_specular_energy(L, mic, src, alpha)
    % Calculates Total Energy (Direct + Reflections) using Image Source Method
    % L: Room Dimensions [Lx Ly Lz]
    % mic: Receiver Pos [x y z]
    % src: Source Pos [x y z]
    % alpha: Wall Absorption Coefficient (Energy)
    
    order = 3; % 3rd Order captures most early energy in dry rooms
    
    % Reflection Coefficient (Energy)
    R_energy = 1 - alpha;
    
    % Direct Path
    d2 = sum((mic - src).^2);
    E_total = 1/d2; 
    
    % Image Sources Loop
    for x = -order:order
        for y = -order:order
            for z = -order:order
                if x==0 && y==0 && z==0, continue; end % Skip direct
                
                % Image Source Position
                % Formula: s_img = [x*2Lx, y*2Ly, z*2Lz] +/- src
                % Correct general ISM logic:
                % Each dimension adds 2*L shift. 
                % Even index: displacement is (2*n*L)
                % Odd index: displacement is (2*n*L) + (L - s) - s ... mirrored
                
                % Simple Box Model Implementation:
                sx = get_img_coord(src(1), L(1), x);
                sy = get_img_coord(src(2), L(2), y);
                sz = get_img_coord(src(3), L(3), z);
                
                dist_sq = (mic(1)-sx)^2 + (mic(2)-sy)^2 + (mic(3)-sz)^2;
                
                % Order of reflection = |x| + |y| + |z|
                N_refl = abs(x) + abs(y) + abs(z);
                
                % Energy Attenuation by walls
                % E = (R^N) / d^2
                E_total = E_total + (R_energy^N_refl) / dist_sq;
            end
        end
    end
end

function coord = get_img_coord(s, L, n)
    % Calculate coordinate of n-th image source
    if mod(n, 2) == 0
        coord = n*L + s;
    else
        coord = n*L + (L - s); % Mirrored
    end
end

% --- FILE HELPER ---
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