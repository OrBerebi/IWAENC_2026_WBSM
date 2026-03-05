function [binaural_sig, binaural_sig_d, binaural_sig_r, coherence_debug] = compass_bsm_encoder_bl(mic_signals, ATF_direct, HRTF_direct, c_BSM, stft_params, source_orig_fs)
%COMPASS_BSM_ENCODER Signal-Dependent BSM (Band-Limited)
%
%   ROBUST VERSION:
%   1. Band-Limits LCMP: Only runs LCMP up to 'source_orig_fs/2'.
%      Above this, it falls back to a stable Matched Filter.
%   2. Prevents "Musical Noise" caused by upsampled empty bands.
%
%   INPUTS:
%     ...
%     source_orig_fs : [Scalar] Original Sampling Rate of source (e.g. 16000).
%                      Default: 48000 (Full Band).
%
%   Reference:
%   [1] Berger et al., "Performance and Robustness of BSM...", EURASIP 2025.

    %% 1. Setup & Defaults
    if nargin < 6, source_orig_fs = 16e3; end % Default to full bandwidth if not specified
    
    [nSamples, Q] = size(mic_signals);
    nfft = stft_params.nfft;
    hop  = stft_params.hop;
    win  = stft_params.win;
    fs_sim = 48000; % Simulation sampling rate
    nBins_target = nfft/2 + 1;
    
    % --- Determine Active Bandwidth ---
    % Calculate the bin index corresponding to the Nyquist of the ORIGINAL signal
    nyquist_orig = source_orig_fs / 2;
    bin_freqs = linspace(0, fs_sim/2, nBins_target);
    
    % Find the cut-off bin (e.g., bin closest to 8000 Hz)
    [~, max_active_bin] = min(abs(bin_freqs - nyquist_orig));
    
    fprintf('COMPASS-BSM: Band-Limiting LCMP to %d Hz (Bin %d/%d)\n', ...
        nyquist_orig, max_active_bin, nBins_target);

    % --- Input Validation & Resampling (Same as before) ---
    resample_data = @(M, nIn, nOut) perform_resampling(M, nIn, nOut);
    
    [Q_v, nBins_v, K] = size(ATF_direct);
    if nBins_v ~= nBins_target
        ATF_direct = resample_data(ATF_direct, nBins_v, nBins_target); 
    end
    
    [nEars, nBins_h, K_h] = size(HRTF_direct);
    if nBins_h ~= nBins_target
        HRTF_direct = resample_data(HRTF_direct, nBins_h, nBins_target); 
    end
    
    [Q_c, nEars_c, nBins_c] = size(c_BSM);
    if nBins_c ~= nBins_target
        % FIX: c_BSM is [Q x 2 x nBins].
        % The resampler expects [d1 x nBins x d3].
        
        % 1. Permute to [Q x nBins x 2]
        c_BSM_temp = permute(c_BSM, [1, 3, 2]);
        
        % 2. Resample
        c_BSM_res = resample_data(c_BSM_temp, nBins_c, nBins_target);
        
        % 3. Permute back to [Q x 2 x nBins]
        c_BSM = permute(c_BSM_res, [1, 3, 2]);
    end
    
    % Fix for swapped ambient ears
    c_BSM = c_BSM(:, [2, 1], :); 

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
            
            % Update Covariance (kept for LCMP)
            cx_inst = x_k * x_k';
            Cx_mem(:, :, k) = alpha * Cx_mem(:, :, k) + (1-alpha) * cx_inst;
            Rx = Cx_mem(:, :, k); 
            
            % -------------------------------------------------------------
            % 1. Source Estimation (Band-Limited LCMP)
            % -------------------------------------------------------------
            
            if k <= max_active_bin
                % --- ACTIVE BAND (<8kHz): Use High-Performance LCMP ---
                tr_Cx = trace(Rx);
                beta = 0.15 * tr_Cx + 1e-7; % Robust Beta
                R_loaded = Rx + beta * eye(Q);
                
                denom = V_d' * (R_loaded \ V_d);
                if rcond(denom) < 1e-6
                     % Fallback if matrix is singular even in active band
                     vec_energies = sum(abs(V_d).^2, 1).';
                     W_d = (1 ./ vec_energies) .* V_d';
                else
                     W_d = denom \ (V_d' / R_loaded); 
                end
            else
                % --- DEAD BAND (>8kHz): Force Stable Matched Filter ---
                % Do NOT use LCMP here. The data is noise.
                % Use Delay-and-Sum (Matched Filter) to just pass whatever 
                % coherent energy might exist (likely none) without blowing up.
                vec_energies = sum(abs(V_d).^2, 1).'; 
                vec_energies(vec_energies < 1e-12) = 1e-6;
                W_d = (1 ./ vec_energies) .* V_d'; 
                
                % Optional: Aggressively attenuate this band if you want strictly clean output
                W_d = W_d * 0.1; 
            end
            
            % Limiter
            if norm(W_d) > 10, W_d = W_d / norm(W_d) * 10; end
            
            % Estimate Sources
            s_est = W_d * x_k; 
            
            % -------------------------------------------------------------
            % 2. Rendering
            % -------------------------------------------------------------
            p_d = h_d * s_est; 
            
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

        
        
        % --- Synthesis ---
        reconstruct = @(BF) real(ifft([BF, conj(BF(:, end-1:-1:2))], nfft, 2));
        
        buf_total(idx_start:idx_end, :)  = buf_total(idx_start:idx_end, :) + reconstruct(Bin_total).';
        buf_direct(idx_start:idx_end, :) = buf_direct(idx_start:idx_end, :) + reconstruct(Bin_dir).';
        buf_reverb(idx_start:idx_end, :) = buf_reverb(idx_start:idx_end, :) + reconstruct(Bin_rev).';
    end

    %analyze_rx_stability(Cx_mem, fs_sim);
     
    binaural_sig   = buf_total(1:nSamples, :);
    binaural_sig_d = buf_direct(1:nSamples, :);
    binaural_sig_r = buf_reverb(1:nSamples, :);
    
    binaural_sig = circshift(binaural_sig,54,1);

    % figure('Name', 'COMPASS-BSM (Band-Limited)', 'Color', 'w');
    % imagesc(squeeze(coherence_debug(:,:,1))); axis xy; colormap jet; colorbar; clim([0 1]);
    % title('Spatial Coherence');
end

function M_out = perform_resampling(M_in, nIn, nOut)
%PERFORM_RESAMPLING_ROBUST Resamples filters via Time-Domain Zero-Padding
%
%   Instead of interpolating frequency bins (which corrupts phase),
%   this method:
%   1. IFFTs the short filter to Time Domain.
%   2. Zero-pads it to the target length (Linear Interpolation in Freq).
%   3. FFTs it back to the new resolution.
%
%   INPUTS:
%     M_in : [Channels x nBins_Old x K] (Complex Frequency Data)
%     nIn  : Old number of bins (e.g., 513)
%     nOut : New number of bins (e.g., 16385)

    if nIn == nOut
        M_out = M_in;
        return;
    end

    % 1. Determine FFT Sizes
    % nBins = N/2 + 1  ->  N = (nBins - 1) * 2
    N_old = (nIn - 1) * 2;
    N_new = (nOut - 1) * 2;
    
    [d1, ~, d3] = size(M_in);
    
    % 2. IFFT to Time Domain (Short)
    % Reconstruct full spectrum for IFFT
    % Input is [Ch x Freq x K]. Permute to [Freq x Ch x K] for IFFT
    M_perm = permute(M_in, [2, 1, 3]); 
    
    % Create symmetric negative frequencies for real-signal IFFT
    M_full = [M_perm; conj(M_perm(end-1:-1:2, :, :))];
    
    % IFFT -> Time Domain Impulse Response [N_old x Ch x K]
    h_time = real(ifft(M_full, N_old, 1));
    
    % 3. Zero-Pad to New Length
    % Just add zeros to the end. This increases frequency resolution.
    h_time_padded = [h_time; zeros(N_new - N_old, d1, d3)];
    
    % 4. FFT to Frequency Domain (Long)
    H_new_full = fft(h_time_padded, N_new, 1);
    
    % 5. Crop to Positive Frequencies (nOut)
    H_new_half = H_new_full(1:nOut, :, :);
    
    % 6. Restore Dimensions [Ch x nBins_New x K]
    M_out = permute(H_new_half, [2, 1, 3]);
    
    % Optional: Renormalize energy if needed (Parseval's theorem), 
    % but typically strictly preserving the IR values is correct for filtering.
end


function analyze_rx_stability(Rx, fs)
%ANALYZE_RX_STABILITY Visualizes the health of the Covariance Matrix
%
%   INPUTS:
%     Rx : [Q x Q x nBins] The spatial covariance matrix
%     fs : Sampling rate (e.g., 48000)

    [Q, ~, nBins] = size(Rx);
    freqs = linspace(0, fs/2, nBins);
    
    % Initialize metrics
    cond_num = zeros(1, nBins);
    eig_vals = zeros(Q, nBins);
    trace_val = zeros(1, nBins);
    
    fprintf('Analyzing Rx Stability across %d bins...\n', nBins);
    
    for k = 1:nBins
        R_k = Rx(:, :, k);
        
        % 1. Condition Number (Stability of Inversion)
        % High number = Unstable / Singular
        cond_num(k) = cond(R_k);
        
        % 2. Eigenvalues (Signal Energy vs Noise Floor)
        e = sort(real(eig(R_k)), 'descend');
        eig_vals(:, k) = e;
        
        % 3. Total Energy
        trace_val(k) = trace(R_k);
    end
    
    % --- PLOTTING ---
    figure('Name', 'Rx Matrix Diagnosis', 'Color', 'w', 'Position', [100 100 1000 800]);
    
    % Subplot 1: Condition Number (The Smoking Gun)
    subplot(3,1,1);
    semilogy(freqs, cond_num, 'LineWidth', 1.5, 'Color', 'r');
    grid on;
    title('Matrix Instability: Condition Number (Log Scale)');
    ylabel('Condition Number (cond)');
    xlabel('Frequency (Hz)');
    xline(8000, '--k', '8kHz (Orig Nyquist)', 'LabelVerticalAlignment', 'bottom');
    xlim([0 fs/2]);
    
    % Interpretation aid
    text(fs/4, max(cond_num)/10, 'Low = Stable', 'Color', 'g', 'FontWeight', 'bold');
    text(fs/4, max(cond_num), 'High = Unstable (Artifacts)', 'Color', 'r', 'FontWeight', 'bold');

    % Subplot 2: Eigenvalue Spectrum
    subplot(3,1,2);
    plot(freqs, 10*log10(eig_vals(1,:)+1e-12), 'k', 'LineWidth', 1.5); hold on;
    for q=2:Q
        plot(freqs, 10*log10(eig_vals(q,:)+1e-12), '--', 'LineWidth', 1);
    end
    grid on;
    title('Eigenvalue Spectrum (Signal vs Noise Floor)');
    ylabel('Power (dB)');
    xlabel('Frequency (Hz)');
    xline(8000, '--k', '8kHz Limit');
    legend('Dominant Eigenvalue (Signal)', 'Noise Eigenvalues');
    xlim([0 fs/2]);
    
    % Subplot 3: Diagonal Loading Requirement
    % Calculates how much "Beta" you need to stabilize the matrix
    subplot(3,1,3);
    target_cond = 100; % We want condition number < 100
    % Approx formula: cond = (lambda_max + beta) / (lambda_min + beta)
    % beta approx (lambda_max - target*lambda_min) / (target - 1)
    lambda_max = eig_vals(1, :);
    lambda_min = eig_vals(end, :);
    needed_beta = (lambda_max - target_cond * lambda_min) / (target_cond - 1);
    needed_beta(needed_beta < 0) = 0; % No loading needed if already stable
    
    plot(freqs, 10*log10(needed_beta ./ (trace_val + 1e-12) + 1e-12), 'b', 'LineWidth', 1.5);
    grid on;
    title('Required Regularization (% of Trace) to Stabilize');
    ylabel('Beta / Trace (dB)');
    xlabel('Frequency (Hz)');
    xline(8000, '--k', '8kHz Limit');
    yline(-20, 'r:', '1% Loading');
    yline(-10, 'm:', '10% Loading');
    xlim([0 fs/2]);
    
end
