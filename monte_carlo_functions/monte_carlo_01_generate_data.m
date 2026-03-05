function monte_carlo_01_generate_data(base_path, drr_val, fov_cond, room_types ,target_counts ,data_dir, scenarios_file)
    % Optimized with Smart Caching & Randomized Sampling.
    

    timit_root = "/Users/orberebi/Documents/work/datasets/converted_TIMIT/";

    % Simulation Parameters
    pct_to_run = 1.0; 
    sim_table_params.N_mc_per_config = 25;
    sim_table_params.base_path = base_path;
    sim_table_params.timit_root = timit_root;
    sim_table_params.room_types = room_types; 
    sim_table_params.target_counts = target_counts;
    
    % Condition Logic
    sim_table_params.fov_condition = fov_cond; % "Inside" or "Outside"
    sim_table_params.snr_range = [18, 20];   % dB
    
    if drr_val == 999 % Infinity DRR
        sim_table_params.diffuese_flag = false;
        sim_table_params.specular_flag = false; 
        sim_table_params.drr_range = [5, 7]; % e.g., for 12dB -> [11, 13]
    else
        sim_table_params.diffuese_flag = true; 
        sim_table_params.specular_flag = true; 
        sim_table_params.drr_range = [drr_val - 1, drr_val + 1]; % e.g., for 12dB -> [11, 13]
    end

    % Initialize Array
    arraysetup.M            = 5;
    arraysetup.arrayType    = 6;
    arraysetup.sphereType   = "rigid";
    [arraysetup.th_array, ~, arraysetup.ph_array, arraysetup.r_array] = ...
        BSM_toolbox.GetArrayPositions(arraysetup.arrayType, arraysetup.M, 0);

    fprintf('\n=== STEP 1: Generating Simulation Table ===\n');
    % IMPORTANT: generate_simulation_table_specular_mc_v2 must save the output to `scenarios_file`
    generate_simulation_table_specular_mc_v2(sim_table_params, scenarios_file); 

    fprintf('\n=== STEP 3: Execution Loop (Save-As-You-Go) ===\n');
    load(scenarios_file, 'scenarioTable');

    if ~exist(data_dir, 'dir'), mkdir(data_dir); end

    N_total = height(scenarioTable);
    N_exec  = ceil(pct_to_run * N_total);
    indices_to_run = sort(randperm(N_total, N_exec));

    babble_dir = fullfile(data_dir, 'babble_base');
    if ~exist(babble_dir, 'dir'), mkdir(babble_dir); end

    total_timer = tic;

    for k = 1:length(indices_to_run)
        loop_timer = tic;
        idx = indices_to_run(k);
        row = scenarioTable(idx, :);
        
        scene_filename = fullfile(data_dir, sprintf('scene_%04d.mat', idx));
        if isfile(scene_filename)
            continue;
        end

        room_type_str = char(row.RoomType);
        roomsetup = setup_room_params(room_type_str, sim_table_params.diffuese_flag, sim_table_params.specular_flag);
        
        scene_data = struct();
        scene_data.ID = row.ID;
        scene_data.Room = row.RoomType;
        
        babble_file = fullfile(babble_dir, sprintf('babble_%s.mat', room_type_str));
        
        if isfile(babble_file)
            cached = load(babble_file, 'base_babble_rir');
            base_babble_rir = cached.base_babble_rir;
        else
            if iscell(row.BabblePos), b_pos = row.BabblePos{1}; else, b_pos = row.BabblePos; end
            roomsetup.sourcePos = b_pos;
            roomsetup.src_gain  = 1.0; 
            base_babble_rir = room_auralisation_v4(roomsetup, arraysetup,'Tom');
            save(babble_file, 'base_babble_rir');
        end
        
        scene_data.Babble_RIR = base_babble_rir;
        scene_data.Babble_RIR.p_array_t = base_babble_rir.p_array_t * row.BabbleGain;
        scene_data.Babble_RIR.RIR_SH_t  = base_babble_rir.RIR_SH_t * row.BabbleGain; 
        scene_data.Babble_RIR.src_gain  = row.BabbleGain; 
        
        clear base_babble_rir cached; 
        
        if iscell(row.TargetPositions), t_pos_mat = row.TargetPositions{1}; else, t_pos_mat = row.TargetPositions; end
        if iscell(row.TargetGains),     t_gains   = row.TargetGains{1};     else, t_gains   = row.TargetGains;     end
        
        scene_data.Target_RIRs = cell(row.NumTargets, 1);
        for t = 1:row.NumTargets
            roomsetup.sourcePos = t_pos_mat(t, :);
            roomsetup.src_gain  = t_gains(t); 
            scene_data.Target_RIRs{t} = room_auralisation_v4(roomsetup, arraysetup,'Tom');
        end
        scene_data.arrayPos = roomsetup.arrayPos;
        scene_data.T60 = roomsetup.T60_target;

        save(scene_filename, 'scene_data', '-v7.3');
        clear scene_data; 
        
        fprintf('Processed Scene %d/%d (Index %d)\n', k, N_exec, idx);
    end
end