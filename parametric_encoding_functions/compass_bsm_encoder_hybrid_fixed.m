function [binaural_sig, binaural_sig_d, binaural_sig_r, coherence_debug] = compass_bsm_encoder_hybrid_fixed(mic_signals, ATF_direct, HRTF_direct, c_BSM, stft_params, source_orig_fs)
% COMPASS_BSM_ENCODER_HYBRID_FIXED
% Hybrid Encoder with Aliasing-Protected BSM Weights.
%
% CHANGES:
% 1. "Golden Method" Weight Preparation:
%    - PCHIP Interpolation.
%    - Time-Domain Global Alignment (Preserves ICPD).
%    - Strict Windowing/Truncation (Prevents STFT Aliasing).
% 2. Hybrid Logic:
%    - Low Freq (< 1.5k): Uses "Cleaned" Static BSM.
%    - High Freq (> 1.5k): Uses Parametric LCMP + "Cleaned" BSM for reverb.

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
    
    % Cutoffs
    nyquist_orig = source_orig_fs / 2;
    [~, max_active_bin] = min(abs(bin_freqs - nyquist_orig));
    
    hybrid_cutoff_freq = 1500; 
    [~, hybrid_idx] = min(abs(bin_freqs - hybrid_cutoff_freq));
    
    fprintf('COMPASS-BSM HYBRID (Fixed):\n');
    fprintf('  - Low Freq (< %.1f Hz): Cleaned Static BSM\n', bin_freqs(hybrid_idx));
    fprintf('  - High Freq (> %.1f Hz): Parametric LCMP\n', bin_freqs(hybrid_idx));

    % --- STANDARD RESAMPLING for Vectors (ATF/HRTF) ---
    % These are short vectors, so standard DFT resampling is usually fine.
    resample_data = @(M, nIn, nOut) perform_resampling(M, nIn, nOut);
    
    [Q_v, nBins_v, K] = size(ATF_direct);
    if nBins_v ~= nBins_target, ATF_direct = resample_data(ATF_direct, nBins_v, nBins_target); end
    
    [nEars, nBins_h, K_h] = size(HRTF_direct);
    if nBins_h ~= nBins_target, HRTF_direct = resample_data(HRTF_direct, nBins_h, nBins_target); end

    %% 2. "GOLDEN METHOD" BSM WEIGHT PREPARATION
    % We process c_BSM specially to prevent Low-Freq Aliasing
    
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
    c_BSM = permute(c_interp, [1, 2, 3]); % [Q x Ears x Freq]
    
    % B. Ear Swap (As requested)
    c_BSM = c_BSM(:, [2, 1], :);

    % C. Time-Limiting & Global Alignment
    filt_len_safe = 256; 
    
    % 1. Get Impulse Responses for ALL channels
    c_temp_perm = permute(c_BSM, [3, 1, 2]); % [Freq x Q x Ears]
    H_full = [c_temp_perm; conj(c_temp_perm(end-1:-1:2, :, :))];
    h_time_all = real(ifft(H_full, nfft, 1)); 
    
    % 2. Find Global Delay (Center of Energy)
    energy_profile = squeeze(sum(sum(abs(h_time_all).^2, 2), 3));
    [~, peak_idx] = max(energy_profile);
    
    % 3. Calculate Shift to Safe Zone (e.g., sample 64)
    target_start = 64; 
    global_shift = target_start - peak_idx;
    
    % 4. Apply Global Shift
    h_shifted = circshift(h_time_all, global_shift, 1);
    
    % 5. Window/Truncate
    window_mask = zeros(nfft, 1);
    window_mask(1:filt_len_safe) = 1;
    h_clean = h_shifted .* window_mask;
    
    % 6. FFT back to Frequency
    H_clean = fft(h_clean, nfft, 1);
    c_BSM = permute(H_clean(1:nBins_target, :, :), [2, 3, 1]);

    %% 3. Initialization
    alpha = 0.9; 
    Cx_mem = zeros(Q, Q, nBins_target);
    
    % Buffers
    padded_len = ceil(nSamples/hop)*hop + nfft;
    mic_signals_pad = [mic_signals; zeros(padded_len - nSamples, Q)];
    
    buf_total  = zeros(padded_len, 2);
    buf_direct = zeros(padded_len, 2);
    buf_reverb = zeros(padded_len, 2);
    
    coherence_debug = zeros(nBins_target, ceil(nSamples/hop), K);
    
    %% 4. Processing Loop
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
            w_bsm = squeeze(c_BSM(:, :, k)); % Using the CLEANED weights
            
            % Covariance Update
            cx_inst = x_k * x_k';
            Cx_mem(:, :, k) = alpha * Cx_mem(:, :, k) + (1-alpha) * cx_inst;
            Rx = Cx_mem(:, :, k); 
            
            % =============================================================
            % HYBRID LOGIC
            % =============================================================
            
            if k < hybrid_idx
                % --- MODE A: STATIC BSM (Low Freq) ---
                % Uses the aliasing-free w_bsm directly.
                p_total = w_bsm' * x_k;
                
                Bin_total(:, k) = p_total;
                Bin_dir(:, k)   = 0; 
                Bin_rev(:, k)   = p_total; 
                
            else
                % --- MODE B: PARAMETRIC (High Freq) ---
                if k <= max_active_bin
                    % LCMP
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
                    % Band-Limit Fallback
                    vec_energies = sum(abs(V_d).^2, 1).'; 
                    vec_energies(vec_energies < 1e-12) = 1e-6;
                    W_d = (1 ./ vec_energies) .* V_d'; 
                    W_d = W_d * 0.1; 
                end
                
                if norm(W_d) > 10, W_d = W_d / norm(W_d) * 10; end
                
                % Estimate Sources
                s_est = W_d * x_k; 
                
                % Rendering
                p_d = h_d * s_est; 
                
                % Reverb (Uses the CLEANED BSM weights)
                x_res = x_k - V_d * s_est;
                p_r = w_bsm' * x_res; 
                
                Bin_total(:, k) = p_d + p_r;
                Bin_dir(:, k)   = p_d;
                Bin_rev(:, k)   = p_r;
                
                % Debug
                if norm(x_k) > 1e-12
                    norm_V = sqrt(sum(abs(V_d).^2, 1)).';
                    coherence_debug(k, frame_idx+1, :) = abs(V_d' * x_k) ./ (norm_V * norm(x_k));
                end
            end
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
    
    % --- Alignment ---
    % Note: The BSM weights were shifted by 'global_shift' to be causal.
    % This means the output (especially low freq) is earlier than before.
    % We apply a heuristic adjustment, but alignment should be verified per case.
    %binaural_sig = circshift(binaural_sig, 54, 1);
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