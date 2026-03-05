function monte_carlo_02_evaluate_data(base_path, data_dir, scenarios_file, results_dir)

    hrtf_path      = fullfile(base_path, 'HRTF', 'HRIR_L2702_koln.sofa');
    babble_dry_path = fullfile(base_path, 'dry_signals', 'babble.wav'); 

    if ~exist(results_dir, 'dir'), mkdir(results_dir); end
    if ~isfile(scenarios_file), error('Scenarios file missing!'); end

    fs_proc     = 16e3;    
    beta        = 0.075;   
    target_lufs = -23;     

    fprintf('\n=== STEP 1: Calculating Spatial Filters ===\n');
    arraysetup.M            = 5;
    arraysetup.arrayType    = 6;
    arraysetup.sphereType   = "rigid";
    arraysetup.fs           = fs_proc; 
    [arraysetup.th_array, ~, arraysetup.ph_array, arraysetup.r_array] = ...
        BSM_toolbox.GetArrayPositions(arraysetup.arrayType, arraysetup.M, 0);

    dummy_params = init_params_local(1, base_path, results_dir);
    dummy_params.fs = fs_proc;
    data_in      = load_and_prepare_data_local(dummy_params); 
    plot_anechoic_analysis = false;
    res_filters  = anechoic_analysis_FOV_BSM_v5(beta,arraysetup, base_path, results_dir, data_in, plot_anechoic_analysis);

    fprintf('\n=== STEP 2: Binaural Rendering & Mixing ===\n');
    load(scenarios_file, 'scenarioTable');
    N_total = height(scenarioTable);

    if isfile(babble_dry_path)
        [babble_dry, fs_b] = audioread(babble_dry_path);
        if fs_b ~= fs_proc, babble_dry = resample(babble_dry, fs_proc, fs_b); end
    else
        babble_dry = pinknoise(10*fs_proc); 
    end

    for i = 1:N_total
        scene_file = fullfile(data_dir, sprintf('scene_%04d.mat', i));
        if ~isfile(scene_file), continue; end
        
        temp = load(scene_file, 'scene_data');
        scene = temp.scene_data;
        row   = scenarioTable(i, :);
        
        fprintf('Processing Scene %d/%d... ', i, N_total);
        
        roomsetup_render = struct();
        roomsetup_render.fs = fs_proc;
        roomsetup_render.HRTF_path = hrtf_path;
        roomsetup_render.sig_path = fullfile(base_path, 'dry_signals', 'casta_short.wav'); 
        
        if iscell(row.AudioFiles), audio_files_wrapper = row.AudioFiles; else, audio_files_wrapper = {row.AudioFiles}; end
        target_paths = audio_files_wrapper{1};
        
        dry_target_sigs = cell(row.NumTargets, 1);
        max_samples = 0;

        DOAs = [];
        ATF_sources = [];
        HRTF_direct = [];
        omega_parametric = [];
        [DOAs.noise.r, DOAs.noise.az, DOAs.noise.el] = calculate_doa(temp.scene_data.arrayPos, row.BabblePos);
        DOAs.noise.idx = find_closest_doa(res_filters.omega_array_WBSM(:,1:2), [DOAs.noise.el,DOAs.noise.az]);
        DOAs.noise.V   = res_filters.V_array_f_WBSM(:,:,DOAs.noise.idx);
        
        for t = 1:row.NumTargets
            t_file = target_paths(t);
            if isfile(t_file)
                [tmp_sig, fs_s] = audioread(t_file);
                if fs_s ~= fs_proc, tmp_sig = resample(tmp_sig, fs_proc, fs_s); end
                dry_target_sigs{t} = tmp_sig;
                if length(tmp_sig) > max_samples, max_samples = length(tmp_sig); end
            else
                dry_target_sigs{t} = [];
            end
        end
        
        if max_samples == 0, max_samples = 5 * fs_proc; end
        len_sig = max_samples;
        
        res_babble_bin = bsm_binaural_renderer_V3(res_filters.BSMobj, res_filters.c_TBSM, res_filters.c_BSM, scene.Babble_RIR, roomsetup_render);
                                              
        if length(babble_dry) < len_sig
            reps = ceil(len_sig / length(babble_dry));
            temp_babble = repmat(babble_dry, reps, 1);
            dry_b = temp_babble(1:len_sig, :);
        else
            dry_b = babble_dry(1:len_sig, :);
        end
        
        noise_mix.ref      = fftfilt(res_babble_bin.BRIR_ref,  dry_b);
        noise_mix.bsm      = fftfilt(res_babble_bin.BRIR_mls_normal,  dry_b);
        noise_mix.wbsm     = fftfilt(res_babble_bin.BRIR_mls_fov, dry_b);
        mic_signals        = fftfilt(scene.Babble_RIR.p_array_t, dry_b);
        mic_signals_speech = zeros(size(mic_signals));
        
        speech_mix = struct('ref', 0, 'bsm', 0, 'wbsm', 0);
        speech_individual(row.NumTargets) = struct('ref', [], 'bsm', [], 'wbsm', []);
        target_brirs(row.NumTargets)      = struct('ref', [], 'bsm', [], 'wbsm', []); 
        
        for t = 1:row.NumTargets
            res_tgt_bin = bsm_binaural_renderer_V3(res_filters.BSMobj, res_filters.c_TBSM, res_filters.c_BSM, scene.Target_RIRs{t}, roomsetup_render);
            target_brirs(t).ref  = res_tgt_bin.BRIR_ref;
            target_brirs(t).bsm  = res_tgt_bin.BRIR_mls_normal;
            target_brirs(t).wbsm = res_tgt_bin.BRIR_mls_fov;
                                               
            dry_s = dry_target_sigs{t};
            if isempty(dry_s)
                dry_s = zeros(len_sig, 1);
            else
                if length(dry_s) < len_sig
                    dry_s = [dry_s; zeros(len_sig - length(dry_s), size(dry_s, 2))];
                else
                    dry_s = dry_s(1:len_sig, :);
                end
            end
            
            sig_ref  = fftfilt(res_tgt_bin.BRIR_ref,  dry_s);
            sig_bsm  = fftfilt(res_tgt_bin.BRIR_mls_normal,  dry_s);
            sig_wbsm = fftfilt(res_tgt_bin.BRIR_mls_fov, dry_s);
            p_array_sig_curr = fftfilt(scene.Target_RIRs{t}.p_array_t, dry_s);

            mic_signals = mic_signals +  p_array_sig_curr;
            mic_signals_speech = mic_signals_speech +  p_array_sig_curr;

            speech_individual(t).ref  = sig_ref;
            speech_individual(t).bsm  = sig_bsm;
            speech_individual(t).wbsm = sig_wbsm;
            
            speech_mix.ref  = speech_mix.ref  + sig_ref;
            speech_mix.bsm  = speech_mix.bsm  + sig_bsm;
            speech_mix.wbsm = speech_mix.wbsm + sig_wbsm;

            [DOAs.target{t}.r, DOAs.target{t}.az, DOAs.target{t}.el] = calculate_doa(temp.scene_data.arrayPos, row.TargetPositions{1}(t,:));
            DOAs.target{t}.idx = find_closest_doa(res_filters.omega_array_WBSM(:,1:2), [DOAs.target{t}.el,DOAs.target{t}.az]);
            DOAs.target{t}.V = res_filters.V_array_f_WBSM(:,:,DOAs.target{t}.idx);
            DOAs.target{t}.H = squeeze(res_filters.H_f(DOAs.target{t}.idx,:,:)).';

            ATF_sources = cat(3,ATF_sources,DOAs.target{t}.V);
            HRTF_direct = cat(3,HRTF_direct,DOAs.target{t}.H);
            omega_parametric = [omega_parametric;[DOAs.target{t}.el,DOAs.target{t}.az]];
        end

        N = 2^10;           
        win_len = floor(N / 2);     
        stft_params.nfft = N;
        active_win = bartlett(win_len); 
        stft_params.win = [active_win; zeros(N - win_len, 1)]; 
        stft_params.hop = floor(win_len / 2); 
        stft_params.fs  = fs_proc;

        final_mix = struct();
        final_mix.ref  = speech_mix.ref  + noise_mix.ref;
        final_mix.bsm  = speech_mix.bsm  + noise_mix.bsm;
        final_mix.wbsm = speech_mix.wbsm + noise_mix.wbsm;
        
        c_BSM = res_filters.c_BSM.f_mls;
        c_BSM = permute(c_BSM,[1,3,2]);

        lcmp_params.J = 8;           
        lcmp_params.beta = 1e-2;     
        lcmp_params.alpha = 0.95;    

        %[final_mix.com, com_d_mix, com_a_mix]        = compass_bsm_encoder_lcmp_adaptive_smoothed_fast(mic_signals, ATF_sources, HRTF_direct, c_BSM, stft_params,lcmp_params);
        [speech_mix.com, com_d_speech, com_a_speech] = compass_bsm_encoder_lcmp_adaptive_smoothed_fast(mic_signals_speech, ATF_sources, HRTF_direct, c_BSM, stft_params,lcmp_params);

        if exist('adjust_loundness', 'file')
            final_mix.ref  = adjust_loundness(final_mix.ref,  fs_proc, target_lufs);
            final_mix.bsm  = adjust_loundness(final_mix.bsm,  fs_proc, target_lufs);
            final_mix.wbsm = adjust_loundness(final_mix.wbsm, fs_proc, target_lufs);
            %final_mix.com  = adjust_loundness(final_mix.com, fs_proc, target_lufs);
            
            speech_mix.ref  = adjust_loundness(speech_mix.ref,  fs_proc, target_lufs);
            speech_mix.bsm  = adjust_loundness(speech_mix.bsm,  fs_proc, target_lufs);
            speech_mix.wbsm = adjust_loundness(speech_mix.wbsm, fs_proc, target_lufs);
            speech_mix.com  = adjust_loundness(speech_mix.com, fs_proc, target_lufs);
        end
        
        save_filename = fullfile(results_dir, sprintf('eval_scene_%04d.mat', i));
        
        eval_data = struct();
        eval_data.ID = i;
        eval_data.Mix = final_mix;
        eval_data.SpeechMixed = speech_mix;
        eval_data.SpeechSeparated = speech_individual;
        eval_data.NoiseOnly  = noise_mix;  
        eval_data.BRIRs.Babble.ref  = res_babble_bin.BRIR_ref;
        eval_data.BRIRs.Babble.bsm  = res_babble_bin.BRIR_mls_normal;
        eval_data.BRIRs.Babble.wbsm = res_babble_bin.BRIR_mls_fov;
        eval_data.BRIRs.Targets = target_brirs;
        eval_data.fs = fs_proc;
        
        save(save_filename, 'eval_data');
        fprintf('Done.\n');
    end
    
    

end

% --- HELPER FUNCTIONS ---
function params = init_params_local(beta, base_path, save_folder_path)
    params.beta   = beta;
    params.base   = base_path;
    params.save   = base_path + save_folder_path + "/beta_" + num2str(beta) + "/";
    params.HRTFpath    = base_path + "/HRTF/HRIR_L2702_koln.sofa";
    params.ATFpath     = base_path + "/ATF/aria_interpolated.mat";
end

function data = load_and_prepare_data_local(params)
    startup_script(params.base)
    hobj       = sofaToEaro(convertStringsToChars(params.HRTFpath));
    omega_hrtf = [hobj.sourceGrid.elevation,hobj.sourceGrid.azimuth,hobj.sourceGrid.r];

    if hobj.fs ~= params.fs
        hobj = hobj.resampleData(params.fs);
        %hobj.fs = params.fs;
    end
    ATF = load(params.ATFpath).aria;
    ATF.IR = ATF.IR_interpolated;
    %ATF.IR = permute(ATF.IR,[1,3,2]);

     if ATF.fs ~= params.fs
        gComDiv = gcd(double(ATF.fs), params.fs);
        p = double(params.fs / gComDiv);
        q = double(ATF.fs / gComDiv);
        newData = zeros(size(ATF.IR,1), ceil(size(ATF.IR,2) * (p/q)), size(ATF.IR,3));
         % Resample the SH RIRs (Optional, but good for completeness)
        for ii=1:size(ATF.IR,1)
           newData(ii,:,:) = resample(squeeze(ATF.IR(ii,:,:)),p,q);
        end
        ATF.IR = newData;

    end
    ATF.IR = permute(ATF.IR,[1,3,2]);
    az_res = 1;

    % Azimuth plane grid
    az_res = 1;
    %phi_az = deg2rad(linspace(0,180,180/az_res).');
    phi_az = deg2rad(linspace(-180,180,360/az_res).');


    th_az  = deg2rad(90*ones(size(phi_az)));
    r_az   = ones(size(phi_az));
    omega_az = [th_az, phi_az, r_az];

    tmp1 = 2^(nextpow2(size(ATF.IR,3)));
    tmp2 = 2^(nextpow2(size(hobj.data,2)));
    nFFT = max([tmp1,tmp2,2^10]);
    %nFFT = 2^10;

    ATF.IR = cat(3, ATF.IR, zeros(size(ATF.IR,1),size(ATF.IR,2),nFFT-size(ATF.IR,3)));
    hobj.data = cat(2,hobj.data,zeros(size(hobj.data,1),nFFT-size(hobj.data,2),size(hobj.data,3)));
    hobj.taps = nFFT;

    data.hobj = hobj;
    data.ATF = ATF;
    data.omega_hrtf = omega_hrtf;
    data.omega_az   = omega_az;
    data.nFFT = nFFT;
    data.fs   = hobj.fs;
end