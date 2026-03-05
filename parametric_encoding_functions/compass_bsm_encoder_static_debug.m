function [binaural_sig, binaural_sig_d, binaural_sig_r] = compass_bsm_encoder_static_debug(mic_signals, ATF_direct, HRTF_direct, c_BSM, stft_params, source_orig_fs)
% COMPASS_BSM_ENCODER_STATIC_FINAL
% A strict replication of the 'bsm_binaural_renderer_V3' physics inside STFT.
%
% THE FIX:
% 1. PCHIP Interpolation (Matches Ref).
% 2. TIME-LIMITING: We IFFT the weights, window them to a short length
%    (e.g., 256 samples), and FFT back. This prevents STFT Aliasing.
% 3. GLOBAL SHIFT: We align the filter peak to the start to preserve causality.

    %% 1. Setup
    if nargin < 6, source_orig_fs = 16e3; end 
    [nSamples, Q] = size(mic_signals);
    nfft = stft_params.nfft; % 4096
    hop  = stft_params.hop;  % 1024 or 2048
    win  = stft_params.win;  % 4096
    fs_sim = 48000; 
    
    nBins_target = nfft/2 + 1;
    filt_len_safe = 256; % The "Short" filter length used in Reference
    
    %% 2. Prepare BSM Weights (The "Golden" Method)
    
    % A. Interpolate in Frequency (PCHIP)
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
    c_BSM = permute(c_interp, [1, 2, 3]); % [Q x Ears x Freq]
    
    % B. Ear Swap (Essential)
    c_BSM = c_BSM(:, [2, 1], :);

    % C. Time-Limiting (The Critical Fix)
    % We convert to time, center the energy, crop the tail, and convert back.
    % This makes the filter physically short enough for STFT.
    
    c_BSM_clean = zeros(size(c_BSM));
    
    % 1. Get Impulse Responses for ALL channels first to find Global Peak
    c_temp_perm = permute(c_BSM, [3, 1, 2]); % [Freq x Q x Ears]
    H_full = [c_temp_perm; conj(c_temp_perm(end-1:-1:2, :, :))];
    h_time_all = real(ifft(H_full, nfft, 1)); 
    
    % 2. Find Global Delay (Preserves Relative Phase)
    energy_profile = squeeze(sum(sum(abs(h_time_all).^2, 2), 3));
    [~, peak_idx] = max(energy_profile);
    
    % 3. Calculate Shift to "Early" Safe Zone (e.g., sample 64)
    % We don't put it at 0 to avoid wrapping the rising slope.
    target_start = 64; 
    global_shift = target_start - peak_idx;
    
    % 4. Apply Shift and Window
    h_shifted = circshift(h_time_all, global_shift, 1);
    
    % Create a rectangular window for the filter
    % Keeps samples 1 to filt_len_safe, zeros the rest.
    window_mask = zeros(nfft, 1);
    window_mask(1:filt_len_safe) = 1;
    
    h_clean = h_shifted .* window_mask;
    
    % 5. FFT back to Frequency
    H_clean = fft(h_clean, nfft, 1);
    c_BSM_clean = permute(H_clean(1:nBins_target, :, :), [2, 3, 1]);
    
    %% 3. Processing Loop (Standard STFT)
    
    % Buffers
    padded_len = ceil(nSamples/hop)*hop + nfft;
    mic_signals_pad = [mic_signals; zeros(padded_len - nSamples, Q)];
    buf_total  = zeros(padded_len, 2);
    
    nFrames = floor((size(mic_signals_pad,1) - nfft) / hop) + 1;
    reconstruct = @(BF) real(ifft([BF, conj(BF(:, end-1:-1:2))], nfft, 2));

    for frame_idx = 0:nFrames-1
        idx_start = frame_idx * hop + 1;
        idx_end = idx_start + nfft - 1;
        
        % Windowing
        x_win = mic_signals_pad(idx_start:idx_end, :) .* win; 
        
        % FFT
        X_f = fft(x_win, nfft, 1); 
        X_half = X_f(1:nBins_target, :).'; 
        
        Bin_total = zeros(2, nBins_target);
        
        for k = 1:nBins_target
            x_k = X_half(:, k);                     
            w_bsm = squeeze(c_BSM_clean(:, :, k)); % Use CLEAN weights
            
            % Apply Weights
            p_total = w_bsm' * x_k;
            Bin_total(:, k) = p_total;
        end

        y_time = reconstruct(Bin_total).'; 
        buf_total(idx_start:idx_end, :) = buf_total(idx_start:idx_end, :) + y_time;
    end

    binaural_sig = buf_total(1:nSamples, :);
    binaural_sig_d = zeros(size(binaural_sig));
    binaural_sig_r = zeros(size(binaural_sig));
    
    % --- ALIGNMENT ---
    % Since we shifted the filters by 'global_shift', we must undo/account for it.
    % Original Reference usually has a bulk delay.
    % Try aligning this result with the reference using xcorr or finddelay.
    % For now, apply the heuristic shift compensated by our manipulation.
    
    % heuristic_shift = 54 (User provided) + global_shift (Our manipulation)
    % Note: Since we moved energy TO the start (negative shift usually), 
    % the output is earlier.
    
    % Actually, just return the raw signal. You should align it externally
    % using estimate_delay_diff(ref, test) for perfect results.
    % binaural_sig = circshift(binaural_sig, 54, 1); 
end