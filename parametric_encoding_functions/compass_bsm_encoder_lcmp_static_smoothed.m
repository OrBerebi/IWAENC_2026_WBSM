function [binaural_sig, binaural_sig_d, binaural_sig_r] = compass_bsm_encoder_lcmp_static_smoothed(mic_signals, ATF_direct, HRTF_direct, c_BSM, stft_params, lcmp_params )
% COMPASS_BSM_ENCODER_LCMP_FULLBAND
%
% IMPLEMENTATION OF "COMPASS" (Parametric BSM) from the Paper.
%
% KEY FEATURES:
% 1. Eq. 33: Rx Estimation via Time-Averaging + Frequency Smoothing.
% 2. Eq. 16: LCMP Weight Calculation (Data-Dependent).
% 3. Eq. 11: Source Estimation using LCMP weights.
% 4. ROBUST RESAMPLING: Fixes for Split Impulse (ATF) and Transient (HRTF).

    %% 1. Setup
    [nSamples, Q] = size(mic_signals);
    nfft = stft_params.nfft;
    hop  = stft_params.hop;
    win  = stft_params.win;
    nBins_target = nfft/2 + 1;
    fs_sim = stft_params.fs;
    
    ATF_direct = helper_resample_split_robust(ATF_direct([6 5 4 3 2 1], :, :), size(ATF_direct,2), nBins_target); 
    HRTF_direct = helper_resample_pad_robust(HRTF_direct, size(HRTF_direct,2), nBins_target); 
    c_BSM = preprocess_bsm_weights(c_BSM, nBins_target, nfft, fs_sim);

    [Q_v, nBins_v, K] = size(ATF_direct);
    [nEars, nBins_h, K_h] = size(HRTF_direct);
    [Q_c, nEars_c, nBins_c] = size(c_BSM);

    %% 2. DATA-DEPENDENT BEAMFORMER DESIGN (Eq. 33 & Eq. 16)
    %fprintf('COMPASS-BSM: Estimating Rx (Eq. 33) and Computing LCMP Weights (Eq. 16)...\n');
    
    % --- Step A: Estimate Rx (Eq. 33) ---
    % 1. Time Averaging
    Rx_global = compute_global_Rx(mic_signals, nfft, hop, win, nBins_target);
    
    % 2. Frequency Smoothing (Eq. 33 "Sum j=-J to J")
    % J controls bandwidth. We use J=6 (6*2+1=13 bins total) as a safe default.
    J_smooth = lcmp_params.J; 
    beta_val = lcmp_params.beta;
    Rx_smoothed = smooth_Rx_frequency(Rx_global, J_smooth);
    
    % --- Step B: Compute LCMP Weights (Eq. 16) ---
    % W_d = (V' * Rx^-1 * V)^-1 * V' * Rx^-1
    W_LCMP = zeros(Q, K, nBins_target);
    
    for k = 1:nBins_target
        Rx_k = squeeze(Rx_smoothed(:, :, k));
        V_d  = squeeze(ATF_direct(:, k, :)); % [Q x K]
        
        % Diagonal Loading for Robust Inversion (Regularization)
        % (Paper assumes white noise, so we add epsilon to diagonal)
        trace_Rx = trace(Rx_k);
        loading = beta_val * trace_Rx + 1e-9;
        Rx_inv = inv(Rx_k + loading * eye(Q));
        
        % Numerator: V' * Rx^-1
        Num = V_d' * Rx_inv; 
        
        % Denominator: (V' * Rx^-1 * V)
        Denom = Num * V_d;
        
        % Solve: Denom \ Num
        % Result is [K x Q] -> Transpose to [Q x K] for storage
        W_k = (Denom \ Num)'; 
        
        W_LCMP(:, :, k) = W_k;
    end

    %% 3. Processing Loop (Rendering)
    % Now we just apply the pre-computed LCMP weights.
    
    padded_len = ceil(nSamples/hop)*hop + nfft;
    mic_signals_pad = [mic_signals; zeros(padded_len - nSamples, Q)];
    
    buf_total  = zeros(padded_len, 2);
    buf_direct = zeros(padded_len, 2);
    buf_reverb = zeros(padded_len, 2);
    coherence_debug = zeros(nBins_target, ceil(nSamples/hop), 1);
    
    nFrames = floor((size(mic_signals_pad,1) - nfft) / hop) + 1;
    
    for frame_idx = 0:nFrames-1
        idx_start = frame_idx * hop + 1;
        idx_end = idx_start + nfft - 1;
        
        x_win = mic_signals_pad(idx_start:idx_end, :) .* win; 
        X_f = fft(x_win, nfft, 1); 
        X_half = X_f(1:nBins_target, :).'; 
        
        Bin_total = zeros(2, nBins_target);
        Bin_dir   = zeros(2, nBins_target);
        Bin_rev   = zeros(2, nBins_target);
        
        for k = 1:nBins_target
            x_k = X_half(:, k);                     
            V_d = squeeze(ATF_direct(:, k, :));     
            h_d = squeeze(HRTF_direct(:, k, :));    
            w_bsm = squeeze(c_BSM(:, :, k));        
            
            % --- APPLY LCMP (Eq. 11) ---
            % Retrieve pre-computed weights
            W_d = squeeze(W_LCMP(:, :, k)); % [Q x K]
            
            % Estimate Source (Eq. 11): s_hat = W' * x
            s_est = W_d' * x_k; 
            %s_est = zeros(size(s_est));
            
            % Render Direct
            p_d = h_d * s_est; 
            
            % Residual (Subtract Direct Estimator)
            x_res = x_k - V_d * s_est; 
            p_r = w_bsm' * x_res; 
            
            % Sum
            Bin_total(:, k) = p_d + p_r;
            Bin_dir(:, k)   = p_d;
            Bin_rev(:, k)   = p_r;
            
        end

        reconstruct = @(BF) real(ifft([BF, conj(BF(:, end-1:-1:2))], nfft, 2));
        %buf_total(idx_start:idx_end, :)  = buf_total(idx_start:idx_end, :) + reconstruct(Bin_total).';
        buf_direct(idx_start:idx_end, :) = buf_direct(idx_start:idx_end, :) + reconstruct(Bin_dir).';
        buf_reverb(idx_start:idx_end, :) = buf_reverb(idx_start:idx_end, :) + reconstruct(Bin_rev).';
    end

    %binaural_sig   = buf_total(1:nSamples, :);
    binaural_sig_d = buf_direct(1:nSamples, :);
    binaural_sig_r = buf_reverb(1:nSamples, :);
    binaural_sig_r = circshift(binaural_sig_r,665,1);
    binaural_sig = binaural_sig_d + binaural_sig_r;
end

% --- HELPERS ---

function Rx_global = compute_global_Rx(x, nfft, hop, win, nBins)
    % Eq. 33 (Part 1): Time Averaging
    % We iterate through the signal to save memory, accumulating Rx.
    [nSamps, Q] = size(x);
    nFrames = floor((nSamps - nfft) / hop) + 1;
    Rx_sum = zeros(Q, Q, nBins);
    
    for i = 0:nFrames-1
        idx = i*hop + 1;
        segment = x(idx:idx+nfft-1, :) .* win;
        X_f = fft(segment, nfft);
        X_h = X_f(1:nBins, :).'; % [Q x nBins]
        
        for k = 1:nBins
            vec = X_h(:, k);
            Rx_sum(:, :, k) = Rx_sum(:, :, k) + (vec * vec');
        end
    end
    Rx_global = Rx_sum / nFrames;
end

function Rx_smooth = smooth_Rx_frequency(Rx_in, J)
    % Eq. 33 (Part 2): Frequency Smoothing
    % "J controls the bandwidth". Moving average of 2*J+1 bins.
    [Q, ~, nBins] = size(Rx_in);
    Rx_smooth = zeros(size(Rx_in));
    
    for k = 1:nBins
        % Determine indices (clamp to edges)
        idx_start = max(1, k - J);
        idx_end   = min(nBins, k + J);
        
        % Average across frequency dimension
        Rx_smooth(:, :, k) = mean(Rx_in(:, :, idx_start:idx_end), 3);
    end
end

function M_out = helper_resample_split_robust(M_in, nIn, nOut)
    if nIn == nOut, M_out = M_in; return; end
    N_old = (nIn - 1) * 2; N_new = (nOut - 1) * 2;
    [d1, ~, d3] = size(M_in);
    M_perm = permute(M_in, [2, 1, 3]); 
    M_full = [M_perm; conj(M_perm(end-1:-1:2, :, :))];
    h_time = real(ifft(M_full, N_old, 1));
    h_centered = fftshift(h_time, 1);
    pad_total = N_new - N_old;
    if d3 > 1, pad_pre = zeros(floor(pad_total/2), d1, d3); pad_post = zeros(ceil(pad_total/2), d1, d3);
    else, pad_pre = zeros(floor(pad_total/2), d1); pad_post = zeros(ceil(pad_total/2), d1); end
    h_pad = [pad_pre; h_centered; pad_post];
    h_final = ifftshift(h_pad, 1);
    H_new_full = fft(h_final, N_new, 1);
    M_out = permute(H_new_full(1:nOut, :, :), [2, 1, 3]);
end

function M_out = helper_resample_pad_robust(M_in, nIn, nOut)
    if nIn == nOut, M_out = M_in; return; end
    N_old = (nIn - 1) * 2; N_new = (nOut - 1) * 2;
    [d1, ~, d3] = size(M_in);
    M_perm = permute(M_in, [2, 1, 3]); 
    M_full = [M_perm; conj(M_perm(end-1:-1:2, :, :))];
    h_time = real(ifft(M_full, N_old, 1));
    if d3 > 1, pad_zeros = zeros(N_new - N_old, d1, d3);
    else, pad_zeros = zeros(N_new - N_old, d1); end
    h_pad = [h_time; pad_zeros];
    %h_pad = circshift(h_pad,-(nIn-1)/2,1);
    H_new_full = fft(h_pad, N_new, 1);
    M_out = permute(H_new_full(1:nOut, :, :), [2, 1, 3]);
end

function c_BSM_out = preprocess_bsm_weights(c_BSM, nBins_target, nfft, fs_sim)
    % Specialized BSM processing: PCHIP -> Global Alignment -> Time Limiting
    [Q_c, nEars_c, nBins_c] = size(c_BSM);
    
    % A. PCHIP Interpolation
    f_orig   = linspace(0, fs_sim/2, nBins_c).';
    f_target = linspace(0, fs_sim/2, nBins_target).';
    c_interp = zeros(Q_c, nEars_c, nBins_target);
    for q = 1:Q_c
        for e = 1:nEars_c
            vec = squeeze(c_BSM(q, e, :));
            c_interp(q, e, :) = interp1(f_orig, vec, f_target, 'pchip');
        end
    end
    c_BSM = c_interp(:, [2, 1], :); % B. Ear Swap

    % C. Time-Limiting & Global Alignment
    filt_len_safe = (nBins_c-1)*2; 
    c_temp_perm = permute(c_BSM, [3, 1, 2]); % [Freq x Q x Ears]
    H_full = [c_temp_perm; conj(c_temp_perm(end-1:-1:2, :, :))];
    h_time_all = real(ifft(H_full, nfft, 1)); 
    window_mask = zeros(nfft, 1);
    window_mask(1:nfft/2) = 1;
    h_time_all = h_time_all .* window_mask;
    
    % Find Global Delay
    energy_profile = squeeze(sum(sum(abs(h_time_all).^2, 2), 3));
    [~, peak_idx] = max(energy_profile);
    
    % Shift to Safe Zone
    target_start = round(filt_len_safe/8); 
    global_shift = target_start - peak_idx;
    h_shifted = circshift(h_time_all, global_shift, 1);
    
    % Window/Truncate
    window_mask = zeros(nfft, 1);
    window_mask(1:filt_len_safe) = 1;
    h_clean = h_shifted .* window_mask;
    
    % FFT back
    H_clean = fft(h_clean, nfft, 1);
    c_BSM_out = permute(H_clean(1:nBins_target, :, :), [2, 3, 1]);
end