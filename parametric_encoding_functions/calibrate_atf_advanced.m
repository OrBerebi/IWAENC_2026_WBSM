function calibrated_indices = calibrate_atf_advanced(mic_signals, ATF_direct, stft_params, fs)
% CALIBRATE_ATF_ADVANCED
% Checks for Revered Ordering (CW vs CCW) and Phase Rotation.

    [nSamples, Q] = size(mic_signals);
    nfft = stft_params.nfft;
    win  = stft_params.win;
    nBins = nfft/2 + 1;
    
    % --- 1. Get Data Manifold (Principal Eigenvector at ~1.5 kHz) ---
    target_freq = 1500;
    [~, bin_idx] = min(abs(linspace(0, fs/2, nBins) - target_freq));
    
    Cx = zeros(Q, Q);
    mic_pad = [mic_signals; zeros(nfft, Q)];
    nFrames = min(50, floor((nSamples-nfft)/stft_params.hop)); 
    
    for i = 0:nFrames-1
        idx = i*stft_params.hop + 1;
        x_win = mic_pad(idx:idx+nfft-1, :) .* win;
        X_f = fft(x_win, nfft, 1);
        x_k = X_f(bin_idx, 1:Q).'; 
        Cx = Cx + (x_k * x_k');
    end
    
    [E, ~] = eig(Cx);
    v_data = E(:, end); 
    v_data = v_data / norm(v_data);
    
    % --- 2. Get Model Vector ---
    v_model = squeeze(ATF_direct(:, bin_idx, :)); 
    v_model = v_model / norm(v_model);

    fprintf('Advanced Calibration at %d Hz...\n', target_freq);
    best_score = 0;
    best_perm  = 1:Q;
    best_type  = 'None';
    
    % --- TEST 1: Standard Shifts (CCW) ---
    for s = 0:Q-1
        perm = circshift(1:Q, s);
        v_test = v_model(perm);
        score = abs(v_test' * v_data);
        if score > best_score
            best_score = score;
            best_perm = perm;
            best_type = sprintf('Shift %d (Standard)', s);
        end
    end
    
    % --- TEST 2: Reversed Shifts (CW) ---
    base_rev = Q:-1:1;
    for s = 0:Q-1
        perm = circshift(base_rev, s);
        v_test = v_model(perm);
        score = abs(v_test' * v_data);
        if score > best_score
            best_score = score;
            best_perm = perm;
            best_type = sprintf('Shift %d (Reversed/Clockwise)', s);
        end
    end
    
    fprintf('  Best Match: %s\n', best_type);
    fprintf('  Score: %.3f\n', best_score);
    
    if best_score > 0.9
        fprintf('SUCCESS: Mismatch identified.\n');
        calibrated_indices = best_perm;
    else
        fprintf('FAILURE: Still mismatched. Likely a 90-degree or 30-degree rotation.\n');
        calibrated_indices = 1:Q;
    end
    
    % --- DIAGNOSTIC PLOT: Phase Difference ---
    % If permutations fail, plotting the phase diff reveals the rotation angle.
    figure;
    angle_diff = angle(v_data ./ v_model(best_perm));
    angle_diff = unwrap(angle_diff);
    plot(rad2deg(angle_diff), '-o');
    xlabel('Mic Index'); ylabel('Phase Difference (Deg)');
    title(sprintf('Phase Error Pattern (Best Match: %s)', best_type));
    grid on;
    subtitle('Linear slope = TDOA mismatch. Constant offset = Phase ref mismatch.');
end