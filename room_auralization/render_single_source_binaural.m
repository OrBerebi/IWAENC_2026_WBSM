function render_single_source_binaural(base_path, dry_sig_path, out_dir)
    % A unified script to generate a room auralization for a single source
    % and render four binaural outputs: Ref, BSM, WBSM, and COM.

    close all
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    %% --- 1. Simulation Parameters ---
    beta        = 0.075;   
    target_lufs = -23;
    
    room_type_str = 'medium'; % Example room type
    diffuese_flag = false;
    specular_flag = false;

    roomsetup          = setup_room_params(room_type_str, diffuese_flag, specular_flag);
    fs_proc            = roomsetup.fs; 
    roomsetup.shutup   = false;
    roomsetup.sig_path = dry_sig_path;

    
    hrtf_path = fullfile(base_path, 'HRTF', 'HRIR_L2702_koln.sofa');
    roomsetup.HRTF_path = hrtf_path;
    
    fprintf('\n=== STEP 1: Setting up Array and Room ===\n');
    arraysetup.M            = 5;
    arraysetup.arrayType    = 6;
    arraysetup.sphereType   = "rigid";
    arraysetup.fs           = fs_proc; 
    [arraysetup.th_array, ~, arraysetup.ph_array, arraysetup.r_array] = ...
        BSM_toolbox.GetArrayPositions(arraysetup.arrayType, arraysetup.M, 0);
    
    % Define a single source position [x, y, z] in meters
    th = deg2rad(90);
    ph = deg2rad(45);
    %dis = 0.4;  % DRR 12 dB at medium room 
    dis = 0.85; % DRR 6 dB at medium room 
    %dis = 1.85;  % DRR 0 dB at medium room 
    [src_x,src_y,src_z] = s2c(th,ph,dis); 
    roomsetup.sourcePos = [src_x,src_y,src_z] + roomsetup.arrayPos;
    roomsetup.src_gain  = 1.0;

    %% --- 2. Room Auralisation ---
    fprintf('\n=== STEP 2: Simulating Room Acoustics ===\n');
    % Generate the RIRs for the source
    source_rir = room_auralisation_v4(roomsetup, arraysetup, 'Tom');
    %source_rir = room_auralisation_v4(roomsetup, arraysetup, 'MCRoomSim');
    
    
    % Load and format the dry signal
    [dry_sig, fs_s] = audioread(dry_sig_path);
    if fs_s ~= fs_proc
        dry_sig = resample(dry_sig, fs_proc, fs_s); 
    end
    % Ensure column vector
    if size(dry_sig, 2) > 1, dry_sig = mean(dry_sig, 2); end 


    %% --- 3. Calculating Spatial Filters ---
    fprintf('\n=== STEP 3: Calculating Spatial Filters ===\n');
    dummy_params = init_params_local(beta, base_path, out_dir);
    dummy_params.fs = fs_proc;
    data_in = load_and_prepare_data_local(dummy_params); 
    
    plot_anechoic_analysis = false; % Kept false for standard runs
    res_filters = anechoic_analysis_FOV_BSM_v5(beta, arraysetup, base_path, out_dir, data_in, plot_anechoic_analysis);



    %% --- 4. Binaural Rendering (Ref, BSM, WBSM) ---
    fprintf('\n=== STEP 4: Binaural Rendering & Mixing ===\n');
    
    roomsetup_render = struct();
    roomsetup_render.fs = fs_proc;
    roomsetup_render.HRTF_path = hrtf_path;
    
    % Get BRIRs using the BSM renderer
    res_bin = bsm_binaural_renderer_V3(res_filters.BSMobj, res_filters.c_TBSM, res_filters.c_BSM, source_rir, roomsetup);
    
    % Convolve the dry signal with the rendered BRIRs
    sig_ref  = fftfilt(res_bin.BRIR_ref, dry_sig);
    sig_bsm  = fftfilt(res_bin.BRIR_mls_normal, dry_sig);
    sig_wbsm = fftfilt(res_bin.BRIR_mls_fov, dry_sig);
    %sig_wbsm = fftfilt(res_bin.BRIR_ls_fov, dry_sig);
    
    % Generate microphone signals
    mic_signals = fftfilt(source_rir.p_array_t, dry_sig);

    %% --- 5. COM (Parametric) Rendering ---
    fprintf('\n=== STEP 5: COMPASS Rendering ===\n');
    
    % Calculate DOA for the source relative to the array
    [DOA.r, DOA.az, DOA.el] = calculate_doa(roomsetup.arrayPos, roomsetup.sourcePos);
    DOA.idx = find_closest_doa(res_filters.omega_array_WBSM(:,1:2), [DOA.el, DOA.az]);
    DOA.V = res_filters.V_array_f_WBSM(:,:,DOA.idx);
    DOA.H = squeeze(res_filters.H_f(DOA.idx,:,:)).';
    
    ATF_sources = DOA.V;
    HRTF_direct = DOA.H;
    
    % STFT Parameters
    %N = 2^13;
    N = 2^10;
    win_len = floor(N / 2);     
    stft_params.nfft = N;
    active_win = bartlett(win_len); 
    stft_params.win = [active_win; zeros(N - win_len, 1)]; 
    stft_params.hop = floor(win_len / 2); 
    stft_params.fs  = fs_proc;
    
    c_BSM = permute(res_filters.c_BSM.f_mls, [1, 3, 2]);


    lcmp_params.J = 8;           
    lcmp_params.beta = 1e-2;    
    lcmp_params.alpha = 0.95;  

    [sig_com, sig_com_d, sig_com_r] = compass_bsm_encoder_lcmp_adaptive_smoothed_fast(mic_signals, ATF_sources, HRTF_direct, c_BSM, stft_params, lcmp_params);

    %% --- 6. Final Adjustments and Saving ---
    fprintf('\n=== STEP 6: Formatting and Saving ===\n');
    
    % Package into a final structure
    rendered_audio = struct();
    rendered_audio.ref  = sig_ref;
    rendered_audio.bsm  = sig_bsm;
    rendered_audio.wbsm = sig_wbsm;
    rendered_audio.com  = sig_com;
    
    % Optional loudness normalization
    if exist('adjust_loundness', 'file')
        rendered_audio.ref  = adjust_loundness(rendered_audio.ref,  fs_proc, target_lufs);
        rendered_audio.bsm  = adjust_loundness(rendered_audio.bsm,  fs_proc, target_lufs);
        rendered_audio.wbsm = adjust_loundness(rendered_audio.wbsm, fs_proc, target_lufs);
        rendered_audio.com  = adjust_loundness(rendered_audio.com,  fs_proc, target_lufs);
    end



    % Save as wav
    wav_name = fullfile(out_dir, 'ref.wav');
    audiowrite(wav_name,rendered_audio.ref,fs_proc);
    wav_name = fullfile(out_dir, 'bsm.wav');
    audiowrite(wav_name,rendered_audio.bsm,fs_proc);
    wav_name = fullfile(out_dir, 'wbsm.wav');
    audiowrite(wav_name,rendered_audio.wbsm,fs_proc);
    wav_name = fullfile(out_dir, 'com.wav');
    audiowrite(wav_name,rendered_audio.com,fs_proc);


    % calc metrics
    conf.fs_proc = fs_proc;
    conf.f_low_limit  = 20;
    conf.f_mid_limit  = 1500; 
    conf.f_high_limit = fs_proc/2;

    f_band = [conf.f_mid_limit,conf.f_high_limit];
    ERBstruct = prepare_ERB(f_band, conf.fs_proc);
    ERBstruct.fs = conf.fs_proc;
    ERBstruct.f_band = f_band;

    [err.bsm, ~]  = calculate_all_metrics_stft(rendered_audio.ref, rendered_audio.bsm, fs_proc, conf, ERBstruct);
    [err.wbsm, ~] = calculate_all_metrics_stft(rendered_audio.ref, rendered_audio.wbsm, fs_proc, conf, ERBstruct);
    [err.com, ~]  = calculate_all_metrics_stft(rendered_audio.ref, rendered_audio.com, fs_proc, conf, ERBstruct);


    [bms_ref, ild_ref, ic_ref, f_erb]    = compute_binaural_metrics(rendered_audio.ref, fs_proc, conf);
    [bms_bsm, ild_bsm, ic_bsm, ~]    = compute_binaural_metrics(rendered_audio.bsm, fs_proc, conf);
    [bms_wbsm, ild_wbsm, ic_wbsm, ~] = compute_binaural_metrics(rendered_audio.wbsm, fs_proc, conf);
    [bms_com, ild_com, ic_com, ~]    = compute_binaural_metrics(rendered_audio.com, fs_proc, conf);

    bin_errs.bsm  = compute_binaural_errors(bms_ref, ild_ref, ic_ref, bms_bsm, ild_bsm, ic_bsm);
    bin_errs.wbsm = compute_binaural_errors(bms_ref, ild_ref, ic_ref, bms_wbsm, ild_wbsm, ic_wbsm);
    bin_errs.com  = compute_binaural_errors(bms_ref, ild_ref, ic_ref, bms_com, ild_com, ic_com);

    % Save results
    save_filename = fullfile(out_dir, 'errors.mat');
    save(save_filename,'err' , 'bin_errs', 'res_bin', 'fs_proc');
    
    fprintf('Done! Audio saved to %s\n', save_filename);


    


end

% --- HELPER FUNCTIONS ---
function params = init_params_local(beta, base_path, save_folder_path)
    params.beta   = beta;
    params.base   = base_path;
    params.save   = fullfile(base_path, save_folder_path, sprintf('beta_%s', num2str(beta)));
    params.HRTFpath = fullfile(base_path, 'HRTF', 'HRIR_L2702_koln.sofa');
    params.ATFpath  = fullfile(base_path, 'ATF', 'aria_interpolated.mat');
end

function data = load_and_prepare_data_local(params)
    startup_script(params.base)
    hobj       = sofaToEaro(convertStringsToChars(params.HRTFpath));
    omega_hrtf = [hobj.sourceGrid.elevation, hobj.sourceGrid.azimuth, hobj.sourceGrid.r];

    if hobj.fs ~= params.fs
        hobj = hobj.resampleData(params.fs);
    end
    ATF = load(params.ATFpath).aria;
    ATF.IR = ATF.IR_interpolated;

    if ATF.fs ~= params.fs
        gComDiv = gcd(double(ATF.fs), params.fs);
        p = double(params.fs / gComDiv);
        q = double(ATF.fs / gComDiv);
        newData = zeros(size(ATF.IR,1), ceil(size(ATF.IR,2) * (p/q)), size(ATF.IR,3));
        for ii=1:size(ATF.IR,1)
           newData(ii,:,:) = resample(squeeze(ATF.IR(ii,:,:)), p, q);
        end
        ATF.IR = newData;
    end
    
    ATF.IR = permute(ATF.IR,[1,3,2]);

    az_res = 1;
    phi_az = deg2rad(linspace(-180, 180, 360/az_res).');
    th_az  = deg2rad(90*ones(size(phi_az)));
    r_az   = ones(size(phi_az));
    omega_az = [th_az, phi_az, r_az];

    tmp1 = 2^(nextpow2(size(ATF.IR,3)));
    tmp2 = 2^(nextpow2(size(hobj.data,2)));
    nFFT = max([tmp1, tmp2, 2^10]);

    ATF.IR = cat(3, ATF.IR, zeros(size(ATF.IR,1), size(ATF.IR,2), nFFT-size(ATF.IR,3)));
    hobj.data = cat(2, hobj.data, zeros(size(hobj.data,1), nFFT-size(hobj.data,2), size(hobj.data,3)));
    hobj.taps = nFFT;

    data.hobj = hobj;
    data.ATF = ATF;
    data.omega_hrtf = omega_hrtf;
    data.omega_az   = omega_az;
    data.nFFT = nFFT;
    data.fs   = hobj.fs;
end



function [out, out_f] = calculate_all_metrics_stft(sig_ref, sig_test, fs, conf, ERBstruct)
    win_size = 512; 
    overlap = win_size / 2; % Recommended: 50% overlap for audio metrics
    nfft = 512;
    window = hann(win_size); % Recommended: Hann window to prevent spectral leakage
    
    num_channels = size(sig_ref, 2);
    
    % Get Freq Axis
    [~, f_axis, ~] = spectrogram(sig_ref(:, 1), window, overlap, nfft, fs);
    idx_nmse = find(f_axis >= conf.f_low_limit & f_axis <= conf.f_mid_limit);
    idx_mag  = find(f_axis > conf.f_mid_limit  & f_axis <= conf.f_high_limit);
    
    out_f.f_nmse = f_axis(idx_nmse);
    out_f.f_mag  = f_axis(idx_mag);
    out_f.f_lsd  = f_axis(idx_mag);
    
    nmse_ears = zeros(1, num_channels); 
    mag_ears = zeros(1, num_channels); 
    lsd_ears = zeros(1, num_channels);
    
    nmse_f_ears = zeros(length(idx_nmse), num_channels);
    mag_f_ears  = zeros(length(idx_mag), num_channels);
    lsd_f_ears  = zeros(length(idx_mag), num_channels);
    
    for ch = 1:num_channels
        [S_ref, ~, ~] = spectrogram(sig_ref(:, ch), window, overlap, nfft, fs);
        [S_test, ~, ~] = spectrogram(sig_test(:, ch), window, overlap, nfft, fs);
        ref_sq = abs(S_ref).^2;
        
        % 1. NMSE STFT (Fixed ratio of means)
        err_sq = abs(S_ref - S_test).^2;
        mse_t = sum(err_sq(idx_nmse, :), 1);
        ref_t = sum(ref_sq(idx_nmse, :), 1);
        act = ref_t > 1e-12;
        if any(act)
            %nmse_ears(ch) = mean(mse_t(act) ./ ref_t(act));
            nmse_f_ears(:, ch) = mean(err_sq(idx_nmse, act), 2) ./ (mean(ref_sq(idx_nmse, act), 2) + 1e-12);
            
            nmse_ears(ch) = squeeze(mean(10 * log10(nmse_f_ears(:, ch)+1e-12),1));
        end
        
        % 2. Mag STFT (Fixed ratio of means)
        err_mag = abs(abs(S_ref) - abs(S_test)).^2;
        mse_mag_t = sum(err_mag(idx_mag, :), 1);
        ref_mag_t = sum(ref_sq(idx_mag, :), 1);
        act_m = ref_mag_t > 1e-12;
        if any(act_m)
            %mag_ears(ch) = mean(mse_mag_t(act_m) ./ ref_mag_t(act_m));
            mag_f_ears(:, ch) = mean(err_mag(idx_mag, act_m), 2) ./ (mean(ref_sq(idx_mag, act_m), 2) + 1e-12);
            mag_ears(ch) = squeeze(mean(10 * log10(mag_f_ears(:, ch)+1e-12),1));
        end
        
        % 3. LSD STFT (Valid)
        P_ref  = abs(S_ref(idx_mag, :)).^2 + 1e-12;
        P_test = abs(S_test(idx_mag, :)).^2 + 1e-12;
        log_diff_sq = (10 * log10(P_test ./ P_ref)).^2;
        %lsd_ears(ch) = mean(sqrt(mean(log_diff_sq(:, act_m), 1)));
        lsd_f_ears(:, ch) = sqrt(mean(log_diff_sq(:, act_m), 2));
        lsd_ears(ch) = squeeze(mean(lsd_f_ears(:, ch),1));
    end
    
    out.nmse = mean(nmse_ears) ;
    out_f.nmse = 10 * log10(mean(nmse_f_ears, 2) + 1e-12);
    
    out.mag = mean(mag_ears);
    out_f.mag = 10 * log10(mean(mag_f_ears, 2) + 1e-12);
    
    out.lsd = mean(lsd_ears);
    out_f.lsd = mean(lsd_f_ears, 2);
    
    % 4. BSD STFT (Fixed active frame division)
    % Note: Keeping this zero-overlap to match original explicit block handling.
    num_frames = floor(size(sig_ref, 1) / win_size);
    pad_len = max(win_size, size(ERBstruct.C, 1)); 
    bsd_frame_scores = zeros(num_frames, 1);
    bsd_f_frames = zeros(length(ERBstruct.f_c), num_frames);
    v_frames = 0;
    
    for fr = 1:num_frames
        idx = (fr-1)*win_size + 1 : fr*win_size;
        frame_ref = sig_ref(idx, :);
        
        if sum(frame_ref(:).^2) > 1e-8
            v_frames = v_frames + 1;
            frame_ref_pad  = [frame_ref; zeros(pad_len - win_size, size(frame_ref,2))];
            frame_test_pad = [sig_test(idx, :); zeros(pad_len - win_size, size(sig_test,2))];
            bsd_val = BSD_ERB_fast(frame_ref_pad, frame_test_pad, ERBstruct);
            
            bsd_frame_scores(v_frames) = mean(bsd_val(:)); % Pack tightly using v_frames index
            bsd_f_frames(:, v_frames)  = mean(bsd_val, 2);
        end
    end
    
    if v_frames > 0
        % Calculate using exact v_frames denominator to prevent perfect-score exclusion
        %out.bsd = sum(bsd_frame_scores(1:v_frames)) / v_frames; 
        out_f.bsd = sum(bsd_f_frames(:, 1:v_frames), 2) / v_frames;
        out.bsd   = mean(out_f.bsd,1);
    else
        out.bsd = 0; 
        out_f.bsd = zeros(length(ERBstruct.f_c), 1);
    end
end


function [BMS_erb, ILD_erb, IC_erb, f_erb] = compute_binaural_metrics(p_t, fs, conf)
    % Extracts BMS, ILD, and IC from a binaural signal [N_samples x 2]
    % using 512 samples with no overlap.
    
    win_size = 512;
    overlap = 0;
    nfft = 512;
    
    % Left and Right channel spectrograms (using Rectangular window for no overlap)
    [S_L, F, ~] = spectrogram(p_t(:, 1), rectwin(win_size), overlap, nfft, fs);
    [S_R, ~, ~] = spectrogram(p_t(:, 2), rectwin(win_size), overlap, nfft, fs);
    
    % Expected values (mean across time frames) -> Eq 33
    c_ll = mean(abs(S_L).^2, 2);
    c_rr = mean(abs(S_R).^2, 2);
    c_lr = mean(S_L .* conj(S_R), 2);
    
    % Compute metrics per frequency bin -> Eq 34, 35, 36
    BMS_f = 10 * log10(c_ll + c_rr + 1e-12);
    ILD_f = 10 * log10((c_ll + 1e-12) ./ (c_rr + 1e-12));
    IC_f  = real(c_lr) ./ (sqrt(c_ll .* c_rr) + 1e-12);
    
    % Average over ERB scale 
    erb_scale = 21.4 * log10(4.37 * F / 1000 + 1);
    
    % Define ERB band edges based on conf limits
    num_erb_bands = 32; 
    f_min_erb = 21.4 * log10(4.37 * conf.f_low_limit / 1000 + 1);
    f_max_erb = 21.4 * log10(4.37 * conf.f_high_limit / 1000 + 1);
    erb_edges = linspace(f_min_erb, f_max_erb, num_erb_bands + 1);
    
    BMS_erb = zeros(num_erb_bands, 1);
    ILD_erb = zeros(num_erb_bands, 1);
    IC_erb  = zeros(num_erb_bands, 1);
    f_erb   = zeros(num_erb_bands, 1); 
    
    for b = 1:num_erb_bands
        % Find STFT bins that fall into the current ERB band
        idx = find(erb_scale >= erb_edges(b) & erb_scale < erb_edges(b+1));
        
        if isempty(idx)
            continue; 
        end
        
        % Average the metric values within this ERB band
        BMS_erb(b) = mean(BMS_f(idx));
        ILD_erb(b) = mean(ILD_f(idx));
        IC_erb(b)  = mean(IC_f(idx));
        f_erb(b)   = (10^((erb_edges(b) + erb_edges(b+1))/2 / 21.4) - 1) / 4.37 * 1000;
    end
end

function metrics = compute_binaural_errors(BMS_ref, ILD_ref, IC_ref, BMS_test, ILD_test, IC_test)
    % Computes RMSE and biases across frequency bands
    valid_idx = BMS_ref ~= 0; 
    N_f = sum(valid_idx);
    
    if N_f == 0
        metrics.bms_rmse = NaN;
        metrics.ild_rmse = NaN;
        metrics.ic_rmse  = NaN;
        metrics.ic_bias  = NaN;
        return;
    end
    
    metrics.bms_rmse = sqrt((1/N_f) * sum((BMS_ref(valid_idx) - BMS_test(valid_idx)).^2));
    metrics.ild_rmse = sqrt((1/N_f) * sum((ILD_ref(valid_idx) - ILD_test(valid_idx)).^2));
    metrics.ic_rmse  = sqrt((1/N_f) * sum((IC_ref(valid_idx) - IC_test(valid_idx)).^2));
    
    % Delta IC Bias: IC_test - IC_ref
    % Positive bias = test is too coherent (lacking envelopment)
    metrics.ic_bias = (1/N_f) * sum(IC_test(valid_idx) - IC_ref(valid_idx));
end
