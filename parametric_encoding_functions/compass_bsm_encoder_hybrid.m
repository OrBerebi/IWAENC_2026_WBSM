function [binaural_sig, binaural_sig_d, binaural_sig_r, coherence_debug] = compass_bsm_encoder_hybrid(mic_signals, ATF_direct, HRTF_direct, c_BSM, stft_params, source_orig_fs)
%COMPASS_BSM_ENCODER_HYBRID Signal-Dependent BSM with Low-Freq Fallback
%
%   HYBRID LOGIC:
%   1. F < 1.5 kHz: Signal-Independent BSM (Static Filtering).
%      - Fixes high NMSE at low freqs caused by beamformer ill-conditioning.
%   2. F > 1.5 kHz: Parametric BSM (LCMP + Restoration).
%   3. F > 8.0 kHz: Matched Filter (Band-Limited logic).
%
%   INPUTS:
%     ... same as before ...

    %% 1. Setup & Defaults
    if nargin < 6, source_orig_fs = 16e3; end 
    
    [nSamples, Q] = size(mic_signals);
    nfft = stft_params.nfft;
    hop  = stft_params.hop;
    win  = stft_params.win;
    fs_sim = 48000; 
    nBins_target = nfft/2 + 1;
    
    % --- FREQUENCY BANDS SETUP ---
    bin_freqs = linspace(0, fs_sim/2, nBins_target);
    
    % A. Band-Limiting Cutoff (e.g., 8kHz)
    nyquist_orig = source_orig_fs / 2;
    [~, max_active_bin] = min(abs(bin_freqs - nyquist_orig));
    
    % B. Hybrid Crossover Cutoff (1.5 kHz)
    hybrid_cutoff_freq = 1500; 
    [~, hybrid_idx] = min(abs(bin_freqs - hybrid_cutoff_freq));
    
    fprintf('COMPASS-BSM HYBRID:\n');
    fprintf('  - Low Freq (< %.1f Hz): Static BSM (Signal Independent)\n', bin_freqs(hybrid_idx));
    fprintf('  - High Freq (> %.1f Hz): Parametric BSM (LCMP)\n', bin_freqs(hybrid_idx));
    fprintf('  - Band Limit (> %.1f Hz): Matched Filter Fallback\n', nyquist_orig);

    % --- Input Validation & Resampling (Standard) ---
    resample_data = @(M, nIn, nOut) perform_resampling(M, nIn, nOut);
    
    [Q_v, nBins_v, K] = size(ATF_direct);
    if nBins_v ~= nBins_target, ATF_direct = resample_data(ATF_direct, nBins_v, nBins_target); end
    
    [nEars, nBins_h, K_h] = size(HRTF_direct);
    if nBins_h ~= nBins_target, HRTF_direct = resample_data(HRTF_direct, nBins_h, nBins_target); end
    
    [Q_c, nEars_c, nBins_c] = size(c_BSM);
    if nBins_c ~= nBins_target
        c_BSM_temp = permute(c_BSM, [1, 3, 2]);
        c_BSM_res = resample_data(c_BSM_temp, nBins_c, nBins_target);
        c_BSM = permute(c_BSM_res, [1, 3, 2]);
    end
    c_BSM = c_BSM(:, [2, 1], :); % Swap ears if needed

    %% 2. Initialization
    alpha = 0.9; 
    Cx_mem = zeros(Q, Q, nBins_target);
    
    % Buffers
    padded_len = ceil(nSamples/hop)*hop + nfft;
    mic_signals_pad = [mic_signals; zeros(padded_len - nSamples, Q)];
    
    buf_total  = zeros(padded_len, 2);
    buf_direct = zeros(padded_len, 2);
    buf_reverb = zeros(padded_len, 2);
    
    coherence_debug = zeros(nBins_target, ceil(nSamples/hop), K);
    
    %% 3. Processing Loop
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
            % Load Bin Data
            x_k = X_half(:, k);                     
            V_d = squeeze(ATF_direct(:, k, :));     
            h_d = squeeze(HRTF_direct(:, k, :));    
            w_bsm = squeeze(c_BSM(:, :, k));        
            
            % Update Covariance (Always tracked for consistency, though unused in Low Band)
            cx_inst = x_k * x_k';
            Cx_mem(:, :, k) = alpha * Cx_mem(:, :, k) + (1-alpha) * cx_inst;
            Rx = Cx_mem(:, :, k); 
            
            % =============================================================
            % HYBRID LOGIC START
            % =============================================================
            
            if k < hybrid_idx
                % ---------------------------------------------------------
                % MODE A: Low Frequency (< 1.5 kHz) -> STATIC BSM
                % ---------------------------------------------------------
                % Beamforming is ill-conditioned here. Use standard BSM.
                % y = w_bsm' * x
                
                % Total Signal
                p_total = w_bsm' * x_k;
                
                % Outputs
                Bin_total(:, k) = p_total;
                
                % Note: Signal-Independent BSM doesn't separate Direct/Reverb.
                % We assign all energy to 'Reverb' buffer to represent "Ambience",
                % or you could split it 0/0. Here we assign to Reverb to keep Total consistent.
                Bin_dir(:, k)   = 0; 
                Bin_rev(:, k)   = p_total; 
                
            else
                % ---------------------------------------------------------
                % MODE B: High Frequency (> 1.5 kHz) -> PARAMETRIC
                % ---------------------------------------------------------
                
                % --- 1. Source Estimation ---
                if k <= max_active_bin
                    % Active Band (1.5k - 8k): LCMP
                    tr_Cx = trace(Rx);
                    beta = 0.15 * tr_Cx + 1e-7; 
                    R_loaded = Rx + beta * eye(Q);
                    
                    denom = V_d' * (R_loaded \ V_d);
                    if rcond(denom) < 1e-6
                         vec_energies = sum(abs(V_d).^2, 1).';
                         W_d = (1 ./ vec_energies) .* V_d';
                    else
                         W_d = denom \ (V_d' / R_loaded); 
                    end
                else
                    % Dead Band (> 8k): Matched Filter
                    vec_energies = sum(abs(V_d).^2, 1).'; 
                    vec_energies(vec_energies < 1e-12) = 1e-6;
                    W_d = (1 ./ vec_energies) .* V_d'; 
                    W_d = W_d * 0.1; % Attenuation
                end
                
                if norm(W_d) > 10, W_d = W_d / norm(W_d) * 10; end
                
                % Estimate Sources
                s_est = W_d * x_k; 
                
                % --- 2. Rendering ---
                p_d = h_d * s_est; 
                
                x_res = x_k - V_d * s_est;
                p_r = w_bsm' * x_res; 
                
                Bin_total(:, k) = p_d + p_r;
                Bin_dir(:, k)   = p_d;
                Bin_rev(:, k)   = p_r;
                
                % Debug (Coherence)
                if norm(x_k) > 1e-12
                    norm_V = sqrt(sum(abs(V_d).^2, 1)).';
                    coherence_debug(k, frame_idx+1, :) = abs(V_d' * x_k) ./ (norm_V * norm(x_k));
                end
            end
            % =============================================================
            % HYBRID LOGIC END
            % =============================================================
        end

        % --- Synthesis ---
        reconstruct = @(BF) real(ifft([BF, conj(BF(:, end-1:-1:2))], nfft, 2));
        
        buf_total(idx_start:idx_end, :)  = buf_total(idx_start:idx_end, :) + reconstruct(Bin_total).';
        buf_direct(idx_start:idx_end, :) = buf_direct(idx_start:idx_end, :) + reconstruct(Bin_dir).';
        buf_reverb(idx_start:idx_end, :) = buf_reverb(idx_start:idx_end, :) + reconstruct(Bin_rev).';
    end

    binaural_sig   = buf_total(1:nSamples, :);
    binaural_sig_d = buf_direct(1:nSamples, :);
    binaural_sig_r = buf_reverb(1:nSamples, :);
    
    binaural_sig = circshift(binaural_sig, 54, 1);
end

% --- LOCAL HELPERS ---
function M_out = perform_resampling(M_in, nIn, nOut)
    if nIn == nOut, M_out = M_in; return; end
    N_old = (nIn - 1) * 2;
    N_new = (nOut - 1) * 2;
    [d1, ~, d3] = size(M_in);
    M_perm = permute(M_in, [2, 1, 3]); 
    M_full = [M_perm; conj(M_perm(end-1:-1:2, :, :))];
    h_time = real(ifft(M_full, N_old, 1));
    h_time_padded = [h_time; zeros(N_new - N_old, d1, d3)];
    H_new_full = fft(h_time_padded, N_new, 1);
    H_new_half = H_new_full(1:nOut, :, :);
    M_out = permute(H_new_half, [2, 1, 3]);
end