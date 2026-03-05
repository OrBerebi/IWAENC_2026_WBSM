function diagnose_beamformer_values(mic_signals, ATF_direct, HRTF_direct, stft_params, fs)
% DIAGNOSE_BEAMFORMER_VALUES
% Performs a forensic check on the Matched Filter logic for a SINGLE
% high-energy Time-Frequency bin.
%
% GOAL: Find why s_est is zero.

    fprintf('=== BEAMFORMER SANITY CHECK ===\n');

    % --- 1. Apply Geometry Fix (Standard Protocol) ---
    % reverse_idx = [6, 5, 4, 3, 2, 1];
    % ATF_direct = ATF_direct(reverse_idx, :);
    % fprintf('1. Geometry Correction [6 5 4 3 2 1] applied.\n');

    % --- 2. Setup Dimensions ---
    [nSamples, Q] = size(mic_signals);
    nfft = stft_params.nfft;
    nBins_target = nfft/2 + 1;
    
    % --- 3. PCHIP Resampling (Manual, Safe) ---
    fprintf('2. Resampling Vectors (PCHIP)...\n');
    
    % ATF
    [d1, d2] = size(ATF_direct);
    if d2 ~= nBins_target
        f_in = linspace(0, 1, d2);
        f_out = linspace(0, 1, nBins_target);
        ATF_direct = interp1(f_in, ATF_direct.', f_out, 'pchip').'; % [Q x Freq]
    end
    
    % HRTF
    [d1_h, d2_h] = size(HRTF_direct);
    if d2_h ~= nBins_target
        f_in = linspace(0, 1, d2_h);
        f_out = linspace(0, 1, nBins_target);
        HRTF_direct = interp1(f_in, HRTF_direct.', f_out, 'pchip').'; % [Ears x Freq]
    end

    % --- 4. Find the "Loudest" Bin ---
    % We don't want to debug noise. We want to debug the SOURCE.
    % Compute crude STFT of Center Channel
    x_center = mic_signals(:, 1);
    spec = abs(fft(x_center(1:nfft))); 
    spec = spec(1:nBins_target);
    
    % Ignore DC and very low freqs
    spec(1:10) = 0; 
    [max_val, k_test] = max(spec);
    freq_test = (k_test-1) * fs / nfft;
    
    fprintf('3. Target Frequency: %.1f Hz (Bin %d)\n', freq_test, k_test);
    fprintf('   Signal Magnitude: %.6f\n', max_val);
    
    if max_val < 1e-6
        fprintf('   WARNING: Signal is essentially silence. Debugging noise.\n');
    end

    % --- 5. Extract Vectors for this Bin ---
    x_k = zeros(Q, 1);
    % Grab a frame from the middle of the file to ensure signal presence
    start_samp = floor(nSamples/2);
    x_chunk = mic_signals(start_samp : start_samp+nfft-1, :);
    X_chunk = fft(x_chunk .* stft_params.win, nfft);
    x_k = X_chunk(k_test, :).'; % [Q x 1]
    
    V_d = squeeze(ATF_direct(:, k_test)); % [Q x 1]
    h_d = squeeze(HRTF_direct(:, k_test)); % [2 x 1] or [1 x 1]

    % --- 6. The Forensic Audit ---
    fprintf('\n--- VECTOR VALUES ---\n');
    
    % A. Check Input Signal (x_k)
    mag_x = norm(x_k);
    fprintf('A. Mic Signal |x_k|:    %.6f\n', mag_x);
    if mag_x < 1e-9, fprintf('   FAIL: Mic signal is dead.\n'); end
    
    % B. Check Steering Vector (V_d)
    mag_v = norm(V_d);
    fprintf('B. Steering Vec |V_d|:  %.6f\n', mag_v);
    if mag_v < 1e-9, fprintf('   FAIL: ATF is zero/dead.\n'); end
    if isnan(mag_v), fprintf('   FAIL: ATF contains NaNs.\n'); end
    
    % C. Check HRTF (h_d)
    mag_h = norm(h_d);
    fprintf('C. HRTF |h_d|:          %.6f\n', mag_h);
    if mag_h < 1e-9, fprintf('   FAIL: HRTF is zero. Output will be silenced.\n'); end
    
    % D. Check Phase Alignment (Normalized Dot Product)
    if mag_x > 0 && mag_v > 0
        % V' * x checks alignment. 
        % If V matches x, V'x should be large.
        dot_raw = V_d' * x_k;
        coherence = abs(dot_raw) / (mag_v * mag_x);
        fprintf('D. Alignment (0-1):     %.4f\n', coherence);
        
        if coherence < 0.1
            fprintf('   FAIL: Vectors are Orthogonal (Cancellation).\n');
            fprintf('         Mic Phases: %s\n', mat2str(angle(x_k)', 2));
            fprintf('         ATF Phases: %s\n', mat2str(angle(V_d)', 2));
        elseif coherence > 0.8
            fprintf('   PASS: Good Phase Alignment.\n');
        else
            fprintf('   WARN: Mediocre Alignment.\n');
        end
    end
    
    % --- 7. Simulate Matched Filter ---
    fprintf('\n--- MF CALCULATION ---\n');
    
    % Step 1: Normalize V
    v_ref = abs(V_d(1));
    if v_ref < 1e-12, v_ref = 1e-12; end
    V_norm = V_d / v_ref;
    fprintf('E. Ref Mic Mag (V(1)):  %.6f\n', v_ref);
    
    % Step 2: Compute Weights (Simple MF)
    denom = V_norm' * V_norm;
    W_mf = V_norm / denom;
    fprintf('F. Weight Norm |W|:     %.6f\n', norm(W_mf));
    
    % Step 3: Estimate Source
    s_est = W_mf' * x_k;
    fprintf('G. s_est Magnitude:     %.6f\n', abs(s_est));
    
    % Step 4: Render
    p_d = h_d * s_est;
    fprintf('H. binaural_d Mag:      %.6f\n', norm(p_d));
    
    if norm(p_d) < 1e-6
        fprintf('\nCONCLUSION: OUTPUT IS SILENCE.\n');
        if abs(s_est) < 1e-6
            fprintf('  -> Cause: s_est is zero.\n');
            if coherence < 0.1
                fprintf('  -> Root Cause: Misalignment/Orthogonality.\n');
            elseif mag_x < 1e-6
                fprintf('  -> Root Cause: No Input Signal.\n');
            else
                fprintf('  -> Root Cause: Unknown math error.\n');
            end
        else
             fprintf('  -> Cause: s_est is OK, but HRTF (h_d) is zero.\n');
             fprintf('  -> Check HRTF indexing/resampling.\n');
        end
    else
        fprintf('\nCONCLUSION: OUTPUT SEEMS VALID (Non-Zero).\n');
        fprintf('  If you hear nothing, check: Playback scaling or File writing.\n');
    end
end