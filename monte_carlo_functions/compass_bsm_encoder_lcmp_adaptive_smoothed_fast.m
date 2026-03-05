function [binaural_sig, binaural_sig_d, binaural_sig_r] = compass_bsm_encoder_lcmp_adaptive_smoothed_fast(mic_signals, ATF_direct, HRTF_direct, c_BSM, stft_params, lcmp_params)
% COMPASS_BSM_ENCODER_LCMP_ADAPTIVE_SMOOTHED_V2 (OPTIMIZED)
%
% Optimizations applied:
% 1. Vectorized Temporal Adaptation (Step A) using implicit expansion.
% 2. Vectorized Spectral Smoothing (Step B) using 'movmean'.
% 3. Removed costly 'squeeze' and anonymous functions inside the loop.

    %% 1. Setup

    [nSamples, Q] = size(mic_signals);
    nfft = stft_params.nfft;
    hop  = stft_params.hop;
    win  = stft_params.win;
    nBins_target = nfft/2 + 1;
    fs_sim = stft_params.fs;
    % target_lufs = -23;
    % mic_signals  = adjust_loundness(mic_signals,  fs_sim, target_lufs);

    
    % --- GEOMETRY & RESAMPLING (Robust 3D/Multi-source) ---
    if nBins_target ~= size(ATF_direct, 2)
        ATF_direct = helper_resample_split_robust(ATF_direct, size(ATF_direct,2), nBins_target); 
    end

    if nBins_target ~= size(HRTF_direct, 2)
        HRTF_direct = helper_resample_pad_robust(HRTF_direct, size(HRTF_direct,2), nBins_target); 
    end

    if nBins_target ~= size(c_BSM, 3)
        c_BSM = preprocess_bsm_weights(c_BSM, nBins_target, nfft, fs_sim);
    end
    
    %% 2. Adaptive & Smoothing Parameters
    alpha = lcmp_params.alpha;
    J_smooth = lcmp_params.J;    
    beta_val = lcmp_params.beta;
    
    % Initial Rx memory
    Rx_mem = zeros(Q, Q, nBins_target);
    for k = 1:nBins_target, Rx_mem(:,:,k) = 1e-6 * eye(Q); end

    padded_len = ceil(nSamples/hop)*hop + nfft;
    mic_signals_pad = [mic_signals; zeros(padded_len - nSamples, Q)];

    % --- Output Buffers ---
    buf_direct = zeros(padded_len, 2);
    buf_reverb = zeros(padded_len, 2);

    nFrames = floor((size(mic_signals_pad,1) - nfft) / hop) + 1;
    
    % Pre-allocate and pre-compute loop constants for speed
    I_Q = eye(Q);
    idx_diag = 1 : (Q+1) : (Q^2); % Indices for fast trace computation
    idx_conj = nBins_target-1:-1:2; % Pre-calculate conjugate indices for IFFT

    %% 3. Adaptive Processing Loop
    for frame_idx = 0:nFrames-1
        idx_start = frame_idx * hop + 1;
        idx_end = idx_start + nfft - 1;
        
        % FFT
        X_f = fft(mic_signals_pad(idx_start:idx_end, :) .* win, nfft, 1); 
        X_half = X_f(1:nBins_target, :).'; % Size: [Q x nBins_target]
        
        % --- Step A: Temporal Update (Vectorized) ---
        % Implicit expansion efficiently calculates X * X' for all bins simultaneously
        X_3D = reshape(X_half, Q, 1, nBins_target);
        X_cov = X_3D .* conj(permute(X_3D, [2, 1, 3])); % Yields [Q x Q x nBins_target]
        Rx_mem = alpha * Rx_mem + (1 - alpha) * X_cov;
        
        % --- Step B: Frequency Smoothing (Vectorized) ---
        % movmean with 'shrink' exactly mimics the max(1,k-J):min(N,k+J) logic
        Rx_smoothed = movmean(Rx_mem, 2*J_smooth + 1, 3, 'Endpoints', 'shrink');
        
        % Fast extraction of traces across all bins (Vectorized)
        Rx_flat = reshape(Rx_smoothed, Q^2, nBins_target);
        tr_Rx = sum(Rx_flat(idx_diag, :), 1); % Yields [1 x nBins_target]
        
        Bin_dir   = zeros(2, nBins_target);
        Bin_rev   = zeros(2, nBins_target);
        
        % --- Step C & D: LCMP Filtering ---
        for k = 1:nBins_target
            x_k = X_half(:, k);
            
            % 'reshape' is practically free compared to 'squeeze' in a loop
            V_d = reshape(ATF_direct(:, k, :), Q, []);  % [Q x nSources]
            h_d = reshape(HRTF_direct(:, k, :), 2, []); % [2 x nSources]
            w_bsm = c_BSM(:, :, k);                     % [Q x 2]

            % 1. Robust Rx
            Rx_k = Rx_smoothed(:,:,k) + (beta_val * tr_Rx(k) + 1e-9) * I_Q;
            
            % 2. Compute LCMP Weights (Eq. 16)
            inv_Rx_V = Rx_k \ V_d;
            
            % The M x M constraint matrix
            C_mat = V_d' * inv_Rx_V; 
            
            % Add mild diagonal loading to prevent singularity when sources are highly correlated
            nSources = size(V_d, 2);
            C_mat_reg = C_mat + 1e-7 * eye(nSources); 
            
            % Calculate final weights
            W_d = ( C_mat_reg \ (inv_Rx_V') )';
            
            % 3. Source Estimation (Eq. 11)
            s_est = W_d' * x_k; 
            
            p_d = h_d * s_est;                   % Direct component
            p_r = w_bsm' * (x_k - V_d * s_est);  % Residual/Reverberant component
            
            Bin_dir(:, k)   = p_d;
            Bin_rev(:, k)   = p_r;
        end

        % --- Step E: Reconstruction ---
        % Inlined IFFT reconstruction to avoid anonymous function overhead
        buf_direct(idx_start:idx_end, :) = buf_direct(idx_start:idx_end, :) + real(ifft([Bin_dir, conj(Bin_dir(:, idx_conj))], nfft, 2)).';
        buf_reverb(idx_start:idx_end, :) = buf_reverb(idx_start:idx_end, :) + real(ifft([Bin_rev, conj(Bin_rev(:, idx_conj))], nfft, 2)).';
    end

    binaural_sig_d = buf_direct(1:nSamples, :);
    binaural_sig_r = buf_reverb(1:nSamples, :);
    binaural_sig = binaural_sig_d + binaural_sig_r;

end

% --- UTILS (Robust versions) ---
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
    H_new_full = fft(h_pad, N_new, 1);
    M_out = permute(H_new_full(1:nOut, :, :), [2, 1, 3]);
end

function c_BSM_out = preprocess_bsm_weights(c_BSM, nBins_target, nfft, fs_sim)
    [Q_c, nEars_c, nBins_c] = size(c_BSM);
    
    f_orig   = linspace(0, fs_sim/2, nBins_c).';
    f_target = linspace(0, fs_sim/2, nBins_target).';
    c_interp = zeros(Q_c, nEars_c, nBins_target);
    for q = 1:Q_c
        for e = 1:nEars_c
            vec = squeeze(c_BSM(q, e, :));
            c_interp(q, e, :) = interp1(f_orig, vec, f_target, 'pchip');
        end
    end

    filt_len_safe = (nBins_c-1)*2; 
    c_temp_perm = permute(c_BSM, [3, 1, 2]); 
    H_full = [c_temp_perm; conj(c_temp_perm(end-1:-1:2, :, :))];
    h_time_all = real(ifft(H_full, nfft, 1)); 
    
    energy_profile = squeeze(sum(sum(abs(h_time_all).^2, 2), 3));
    [~, peak_idx] = max(energy_profile);
    
    target_start = round(filt_len_safe/8); 
    global_shift = target_start - peak_idx;
    h_shifted = circshift(h_time_all, global_shift, 1);
    
    window_mask = zeros(nfft, 1);
    window_mask(1:filt_len_safe) = 1;
    h_clean = h_shifted .* window_mask;
    
    H_clean = fft(h_clean, nfft, 1);
    c_BSM_out = permute(H_clean(1:nBins_target, :, :), [2, 3, 1]);
end