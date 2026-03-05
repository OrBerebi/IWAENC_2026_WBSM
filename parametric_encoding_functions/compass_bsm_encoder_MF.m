function [binaural_sig, binaural_sig_d, binaural_sig_r] = compass_bsm_encoder_MF(mic_signals, ATF_direct, HRTF_direct, c_BSM, stft_params)
% COMPASS_BSM_ENCODER_GOLD
%
% FIXES:
% 1. GEOMETRY: ATF Reversed [6 5 4 3 2 1].
% 2. RESAMPLING: Supports K>1 Sources (Fixed 3D Padding).
% 3. BEAMFORMER: Full-Band Matched Filter (Normalized).

    %% 1. Setup
    [nSamples, Q] = size(mic_signals);
    nfft = stft_params.nfft;
    hop  = stft_params.hop;
    win  = stft_params.win;
    fs_sim = stft_params.fs;
    nBins_target = nfft/2 + 1;
    
    % --- GEOMETRY & RESAMPLING (Robust 3D/Multi-source) ---
    ATF_direct  = helper_resample_split_robust(ATF_direct([6 5 4 3 2 1], :, :), size(ATF_direct,2), nBins_target); 
    HRTF_direct = helper_resample_pad_robust(HRTF_direct, size(HRTF_direct,2), nBins_target); 
    c_BSM       = preprocess_bsm_weights(c_BSM, nBins_target, nfft, fs_sim);
    

    %% 2. Processing Loop
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
            V_d = squeeze(ATF_direct(:, k, :));     % [Q x K]
            h_d = squeeze(HRTF_direct(:, k, :));    % [2 x K]
            w_bsm = squeeze(c_BSM(:, :, k));        % [Q x 2]
            
            % --- MATCHED FILTER ---
            

            % (Assuming inputs are already good)
            % 1. Normalize V_d (Corrected for K > 1)
            % We normalize each source column by its own reference mic magnitude.
            v_mag = abs(V_d(1, :)); % [1 x K] vector
            v_mag(v_mag < 1e-12) = 1e-12; % Safety
            
            % FIX: Use Element-wise division (./) with implicit expansion
            %V_steer = V_d ./ v_mag;

            V_steer = V_d;
            
            % Weights: [Q x K]
            denom = V_steer' * V_steer;
            W_d = V_steer / (denom + 1e-12 * eye(size(denom))); % Regularize for K>1
            
            % Estimate: [K x 1]
            s_est = W_d' * x_k; 
            %s_est = zeros(size(s_est));
            
            % Render: [2 x K] * [K x 1] = [2 x 1]
            p_d = h_d * s_est; 
            
            
            % Residual
            x_res = x_k - V_steer * s_est; 
            p_r = w_bsm' * x_res; 
            
            % Sum (Both are 2x1)
            Bin_total(:, k) = p_d + p_r;
            Bin_dir(:, k)   = p_d;
            Bin_rev(:, k)   = p_r;
            
            if norm(x_k) > 1e-12
                % Debug (Just Source 1)
                v1 = V_d(:, 1);
                norm_V = norm(v1);
                coherence_debug(k, frame_idx+1, :) = abs(v1' * x_k) ./ (norm_V * norm(x_k));
            end
        end

        reconstruct = @(BF) real(ifft([BF, conj(BF(:, end-1:-1:2))], nfft, 2));
        %buf_total(idx_start:idx_end, :)  = buf_total(idx_start:idx_end, :) + reconstruct(Bin_total).';
        buf_direct(idx_start:idx_end, :) = buf_direct(idx_start:idx_end, :) + reconstruct(Bin_dir).';
        buf_reverb(idx_start:idx_end, :) = buf_reverb(idx_start:idx_end, :) + reconstruct(Bin_rev).';
    end

    
    binaural_sig_d = buf_direct(1:nSamples, :);
    binaural_sig_r = buf_reverb(1:nSamples, :);
    binaural_sig_r = circshift(binaural_sig_r,665,1);
    binaural_sig = binaural_sig_d + binaural_sig_r;
    
end

% --- ROBUST HELPERS (K > 1 Compatible) ---


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