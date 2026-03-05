function diagnose_beamformer_full_spectrum(mic_signals, ATF_direct, HRTF_direct, stft_params, fs)
% DIAGNOSE_BEAMFORMER_FULL_SPECTRUM_V4
%
% VERIFICATION OF CONTEXT-AWARE RESAMPLING:
% 1. Geometry: Fixed ([6 5 4 3 2 1]).
% 2. ATF Resampling: Split-Impulse Logic (Preserves t=0 wrap).
% 3. HRTF Resampling: Simple-Padding Logic (Preserves t=512 latency).

    fprintf('=== FULL SPECTRUM DIAGNOSIS (V4 - Context-Aware Fix) ===\n');

    % --- 1. Geometry Correction ---
    % reverse_idx = [6, 5, 4, 3, 2, 1];
    % ATF_direct = ATF_direct(reverse_idx, :);
    % fprintf('1. Geometry Correction [6 5 4 3 2 1] applied.\n');

    % --- 2. Setup & RESAMPLING (The Context-Aware Fix) ---
    [nSamples, Q] = size(mic_signals);
    nfft = stft_params.nfft;
    nBins = nfft/2 + 1;
    freqs = linspace(0, fs/2, nBins);
    
    fprintf('2. Resampling Vectors (Context-Aware)...\n');
    
    % A. ATF: Split Impulse Logic
    resample_split = @(M, nIn, nOut) helper_resample_split(M, nIn, nOut);
    if size(ATF_direct, 2) ~= nBins
        ATF_direct = resample_split(ATF_direct, size(ATF_direct, 2), nBins);
    end
    
    % B. HRTF: Simple Padding Logic
    resample_pad = @(M, nIn, nOut) helper_resample_pad(M, nIn, nOut);
    if size(HRTF_direct, 2) ~= nBins
        HRTF_direct = resample_pad(HRTF_direct, size(HRTF_direct, 2), nBins);
    end

    % --- 3. Compute Data Covariance ---
    fprintf('3. Accumulating Covariance stats...\n');
    Cx_sum = zeros(Q, Q, nBins);
    n_frames = floor((nSamples - nfft) / stft_params.hop);
    win = stft_params.win;
    count = 0;
    scan_frames = min(50, n_frames); 
    
    for i = 0:scan_frames-1
        idx = i * stft_params.hop + 1;
        x_win = mic_signals(idx:idx+nfft-1, :) .* win;
        X_f = fft(x_win, nfft);
        for k = 1:nBins
            x_k = X_f(k, :).';
            Cx_sum(:, :, k) = Cx_sum(:, :, k) + (x_k * x_k');
        end
        count = count + 1;
    end
    Cx_sum = Cx_sum / count;

    % --- 4. Spectral Analysis Loop ---
    fprintf('4. Analyzing Bins...\n');
    
    metric_align = zeros(1, nBins);
    pow_input    = zeros(1, nBins); 
    pow_est      = zeros(1, nBins); 
    pow_output   = zeros(1, nBins); 
    norm_atf     = zeros(1, nBins);
    norm_hrtf    = zeros(1, nBins);
    
    for k = 2:nBins % Skip DC
        Rx = Cx_sum(:, :, k);
        V_d = squeeze(ATF_direct(:, k));
        h_d = squeeze(HRTF_direct(:, k));
        
        % A. Alignment Score
        [E, ~] = eig(Rx);
        v_data = E(:, end); 
        num = abs(V_d' * v_data);
        den = norm(V_d) * norm(v_data);
        metric_align(k) = num / (den + 1e-12);
        
        % B. Matched Filter
        % Normalize V to Ref Mic Mag
        v_ref = abs(V_d(1)); if v_ref < 1e-12, v_ref = 1e-12; end
        V_steer = V_d / v_ref;
        
        w_denom = V_steer' * V_steer;
        W = V_steer / (w_denom + 1e-12);
        
        % C. Power Estimates
        pow_input(k)  = real(trace(Rx) / Q);
        pow_est(k)    = real(W' * Rx * W);
        pow_output(k) = (norm(h_d)^2) * pow_est(k);
        
        % D. Norms
        norm_atf(k)  = norm(V_d);
        norm_hrtf(k) = norm(h_d);
    end

    % --- 5. Plotting ---
    figure('Name', 'Full Spectrum Diagnosis V4', 'Color', 'w');
    
    % Plot 1: Alignment
    subplot(3,1,1);
    plot(freqs, metric_align, 'b', 'LineWidth', 1.5);
    yline(0.9, 'g--', 'Excellent');
    yline(0.5, 'r--', 'Poor');
    xlim([0 20000]); ylim([0 1.1]);
    grid on;
    title('1. Geometry Alignment Score (0-1)');
    
    % Plot 2: Signal Levels (dB)
    subplot(3,1,2);
    semilogx(freqs, 10*log10(pow_input + 1e-12), 'k', 'LineWidth', 1, 'DisplayName', 'Input (Avg Mic)');
    hold on;
    semilogx(freqs, 10*log10(pow_est + 1e-12), 'b', 'LineWidth', 1.5, 'DisplayName', 'Est. Source (s_{est})');
    semilogx(freqs, 10*log10(pow_output + 1e-12), 'r', 'LineWidth', 1.5, 'DisplayName', 'Binaural Out (p_d)');
    xlim([100 20000]);
    grid on;
    legend('Location', 'best');
    title('2. Signal Power Levels (dB)');
    subtitle('Fix Check: Red line must follow Blue line (High Freqs match).');
    ylabel('Power (dB)');
    
    % Plot 3: Scaling Factors
    subplot(3,1,3);
    semilogx(freqs, 20*log10(norm_atf + 1e-12), 'b--', 'DisplayName', 'ATF Norm');
    hold on;
    semilogx(freqs, 20*log10(norm_hrtf + 1e-12), 'r--', 'DisplayName', 'HRTF Norm');
    yline(0, 'k-', 'Unit Gain');
    xlim([100 20000]);
    grid on;
    legend;
    title('3. Model Vector Scaling (dB)');
    subtitle('Fix Check: Red line MUST be smooth (No jagged oscillations).');
    ylabel('Magnitude (dB)');
    
    fprintf('Done. Check figures.\n');
end

% --- LOCAL HELPERS ---

function M_out = helper_resample_split(M_in, nIn, nOut)
    % FOR ATF: Preserves wrapping at t=0.
    % Rotate -> Pad -> Unrotate
    if nIn == nOut, M_out = M_in; return; end
    N_old = (nIn - 1) * 2;
    N_new = (nOut - 1) * 2;
    [d1, ~] = size(M_in);
    
    M_perm = permute(M_in, [2, 1]); 
    M_full = [M_perm; conj(M_perm(end-1:-1:2, :))];
    
    h_time = real(ifft(M_full, N_old, 1));
    h_centered = fftshift(h_time, 1);
    
    pad_total = N_new - N_old;
    h_pad = [zeros(floor(pad_total/2), d1); h_centered; zeros(ceil(pad_total/2), d1)];
    
    h_final = ifftshift(h_pad, 1);
    H_new_full = fft(h_final, N_new, 1);
    M_out = permute(H_new_full(1:nOut, :), [2, 1]);
end

function M_out = helper_resample_pad(M_in, nIn, nOut)
    % FOR HRTF: Preserves transient position.
    % Simple Zero Pad at End.
    if nIn == nOut, M_out = M_in; return; end
    N_old = (nIn - 1) * 2;
    N_new = (nOut - 1) * 2;
    [d1, ~] = size(M_in);
    
    M_perm = permute(M_in, [2, 1]); 
    M_full = [M_perm; conj(M_perm(end-1:-1:2, :))];
    
    h_time = real(ifft(M_full, N_old, 1));
    h_pad = [h_time; zeros(N_new - N_old, d1)];
    
    H_new_full = fft(h_pad, N_new, 1);
    M_out = permute(H_new_full(1:nOut, :), [2, 1]);
end