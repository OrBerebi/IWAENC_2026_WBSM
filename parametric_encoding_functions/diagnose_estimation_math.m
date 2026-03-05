function diagnose_estimation_math(mic_signals, ATF_direct, stft_params, fs)
% DIAGNOSE_ESTIMATION_MATH
% Performs mathematical analysis of the Source Estimation problem.
%
% METRICS:
% 1. Spectral Coherence: How well v_model matches v_data shape (0 to 1).
% 2. Magnitude Ratio: |v_model| / |v_data|. Checks for scaling mismatches.
% 3. White Noise Gain: 10*log10(|w|^2). Checks for noise amplification.

    % --- 1. Apply Known Geometry Fix ---
    % We established this is necessary for the vectors to align.
    %reverse_idx = [6, 5, 4, 3, 2, 1];
    %ATF_direct = ATF_direct(reverse_idx, :);

    % --- 2. Setup ---
    [nSamples, Q] = size(mic_signals);
    nfft = stft_params.nfft;
    win  = stft_params.win;
    nBins = nfft/2 + 1;
    freqs = linspace(0, fs/2, nBins);
    
    % Resample ATF if needed
    [~, nBins_atf] = size(ATF_direct);
    if nBins_atf ~= nBins
        % Simple resampling for diagnostic (assumes vectors)
        ATF_perm = permute(ATF_direct, [2, 1]);
        ATF_res = interp1(linspace(0, 1, nBins_atf), ATF_perm, linspace(0, 1, nBins), 'pchip').';
        ATF_direct = ATF_res;
    end

    % --- 3. Extract "Ground Truth" Data Vectors ---
    % We compute the Principal Eigenvector of the Data Covariance
    % averaged over the active portion of the signal.
    
    fprintf('Computing Data Manifold...\n');
    Cx_sum = zeros(Q, Q, nBins);
    
    % Use first 100 frames (approx 50-100ms) where direct sound dominates
    nFrames = min(100, floor((nSamples-nfft)/stft_params.hop)); 
    mic_pad = [mic_signals; zeros(nfft, Q)];
    
    for i = 0:nFrames-1
        idx = i*stft_params.hop + 1;
        x_win = mic_pad(idx:idx+nfft-1, :) .* win;
        X_f = fft(x_win, nfft, 1);
        X_half = X_f(1:nBins, :).';
        
        for k = 1:nBins
            x_k = X_half(:, k);
            Cx_sum(:, :, k) = Cx_sum(:, :, k) + (x_k * x_k');
        end
    end
    
    % --- 4. Compute Metrics per Bin ---
    coherence_score = zeros(1, nBins);
    mag_ratio       = zeros(1, nBins);
    wng_db          = zeros(1, nBins);
    
    for k = 2:nBins % Skip DC
        % A. Get Data Vector (v_d)
        [E, D] = eig(Cx_sum(:, :, k));
        [~, max_idx] = max(diag(D));
        v_data = E(:, max_idx);
        
        % B. Get Model Vector (v_m)
        v_model = squeeze(ATF_direct(:, k));
        
        % Metric 1: Spectral Coherence (Shape Match)
        % |v_m' * v_d| / (|v_m|*|v_d|)
        % 1.0 = Perfect Phase Alignment
        num = abs(v_model' * v_data);
        den = norm(v_model) * norm(v_data);
        coherence_score(k) = num / (den + 1e-12);
        
        % Metric 2: Magnitude Ratio (Scaling Match)
        % |v_model_avg| / |v_data_avg|
        % If this is 0.01, Model is much quieter than Data (Gain Risk)
        mag_model = mean(abs(v_model));
        mag_data  = mean(abs(v_data)); % Eigenvectors are unit norm, so verify scaling
        
        % Actually, let's compare V_model norm to "Unit".
        % Data vector v_data is Unit Norm by definition of eig().
        % So we just look at norm(v_model).
        mag_ratio(k) = norm(v_model); 
        
        % Metric 3: White Noise Gain (WNG)
        % For a Matched Filter: w = v_model / ||v_model||^2
        % WNG = ||w||^2 = 1 / ||v_model||^2
        w_mf = v_model / (norm(v_model)^2 + 1e-12);
        wng_db(k) = 10 * log10(norm(w_mf)^2);
    end
    
    % --- 5. Visualization ---
    figure('Color','w', 'Name', 'Mathematical Estimation Diagnosis');
    
    % Plot 1: Coherence (The "Quality" Check)
    subplot(3,1,1);
    plot(freqs, coherence_score, 'b', 'LineWidth', 1.5);
    yline(0.9, 'g--', 'Good');
    yline(0.5, 'r--', 'Poor');
    xlim([0 fs/2]); ylim([0 1.1]);
    title('1. Spectral Coherence (Direction Match)');
    ylabel('Normalized Dot Product');
    grid on;
    legend('Model vs Data Match');
    
    % Plot 2: Magnitude / WNG (The "Noise" Check)
    subplot(3,1,2);
    plot(freqs, wng_db, 'r', 'LineWidth', 1.5);
    yline(0, 'k--');
    yline(10, 'm--', 'High Gain Risk (>10dB)');
    xlim([0 fs/2]);
    title('2. White Noise Gain (Amplification Risk)');
    ylabel('Gain (dB)');
    grid on;
    legend('Required Gain to Invert ATF');
    
    % Plot 3: Magnitude Spectra (The "Scaling" Check)
    subplot(3,1,3);
    plot(freqs, 20*log10(mag_ratio), 'k', 'LineWidth', 1.5);
    yline(0, 'b--', 'Unit Norm');
    xlim([0 fs/2]);
    title('3. Model Vector Norm (dB)');
    ylabel('Norm (dB re 1.0)');
    grid on;
    xlabel('Frequency (Hz)');
    
    % --- 6. Text Diagnosis ---
    avg_coh = mean(coherence_score(freqs>100 & freqs<8000));
    avg_wng = mean(wng_db(freqs>100 & freqs<8000));
    
    fprintf('=== MATHEMATICAL DIAGNOSIS ===\n');
    fprintf('Avg Coherence (100Hz-8kHz): %.3f  (Should be > 0.9)\n', avg_coh);
    fprintf('Avg WNG (100Hz-8kHz):       %.1f dB (Should be < 10 dB)\n', avg_wng);
    
    if avg_wng > 15
        fprintf('CRITICAL: High White Noise Gain detected.\n');
        fprintf('  The ATF vectors have very small magnitudes (High attenuation).\n');
        fprintf('  Inverting them causes massive noise amplification.\n');
        fprintf('  FIX: Use Unity-Gain (Phase-Only) Beamforming.\n');
    elseif avg_coh < 0.8
        fprintf('CRITICAL: Low Coherence detected.\n');
        fprintf('  Even with geometry fix, vectors do not match well.\n');
        fprintf('  Check for Reflection mismatch or Near-field vs Far-field issues.\n');
    else
        fprintf('PASS: Math looks healthy. Issue might be in the loop implementation.\n');
    end
end