function calibrated_indices = calibrate_atf_geometry(mic_signals, ATF_direct, stft_params, fs)
% CALIBRATE_ATF_GEOMETRY
% Brute-forces Microphone Permutations to maximize alignment score.
%
% Returns the index order that fixes the mismatch (if it's an indexing issue).

    [nSamples, Q] = size(mic_signals);
    nfft = stft_params.nfft;
    win  = stft_params.win;
    nBins = nfft/2 + 1;
    
    % --- 1. Get Data Manifold (Principal Eigenvector at ~2kHz) ---
    % We pick a frequency where the mismatch is visible but not nulled (e.g. 2kHz)
    target_freq = 2000;
    [~, bin_idx] = min(abs(linspace(0, fs/2, nBins) - target_freq));
    
    % Estimate Covariance at this bin
    Cx = zeros(Q, Q);
    mic_pad = [mic_signals; zeros(nfft, Q)];
    nFrames = min(50, floor((nSamples-nfft)/stft_params.hop)); % Use onset
    
    for i = 0:nFrames-1
        idx = i*stft_params.hop + 1;
        x_win = mic_pad(idx:idx+nfft-1, :) .* win;
        X_f = fft(x_win, nfft, 1);
        x_k = X_f(bin_idx, 1:Q).'; % [Q x 1]
        Cx = Cx + (x_k * x_k');
    end
    
    [E, ~] = eig(Cx);
    v_data = E(:, end); % Dominant eigenvector (Data)
    v_data = v_data / norm(v_data);
    
    % --- 2. Get Model Vector ---
    % Ensure ATF is correct shape
    [d1, d2, d3] = size(ATF_direct);
    if d2 ~= nBins % Resample if needed
        % (Simplified resampling for just one bin)
        % Assuming ATF is [Q x nBins] here for simplicity of the check
    end
    v_model = squeeze(ATF_direct(:, bin_idx, :)); 
    v_model = v_model / norm(v_model);

    % --- 3. Test Permutations ---
    fprintf('Testing Microphone Permutations at %d Hz...\n', target_freq);
    
    base_score = abs(v_model' * v_data);
    fprintf('  Baseline Score: %.3f\n', base_score);
    
    best_score = base_score;
    best_perm  = 1:Q;
    
    % A. Circular Shifts (Rotation of circular array)
    for s = 1:Q-1
        perm = circshift(1:Q, s);
        v_shifted = v_model(perm);
        score = abs(v_shifted' * v_data);
        if score > best_score
            best_score = score;
            best_perm = perm;
            fprintf('  Found Better Shift! (%d): Score %.3f\n', s, score);
        end
    end
    
    % B. Reverse Order (Clockwise vs CCW)
    perm_rev = Q:-1:1;
    v_rev = v_model(perm_rev);
    score = abs(v_rev' * v_data);
    if score > best_score
        best_score = score;
        best_perm = perm_rev;
        fprintf('  Found Reverse Match! Score %.3f\n', score);
    end
    
    % Result
    if best_score > 0.95
        fprintf('\nSUCCESS: Found correction! Use these indices: %s\n', mat2str(best_perm));
        calibrated_indices = best_perm;
    else
        fprintf('\nFAILURE: Permutations did not fix it. Likely a Rotation/Angle Mismatch.\n');
        calibrated_indices = 1:Q;
    end
end