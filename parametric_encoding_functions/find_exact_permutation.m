function calibrated_indices = find_exact_permutation(mic_signals, ATF_direct, stft_params, fs)
% FIND_EXACT_PERMUTATION
% Tests all 6 circular shifts of the REVERSED (Clockwise) winding.
% One of these will be the perfect geometric match.

    [nSamples, Q] = size(mic_signals);
    nfft = stft_params.nfft;
    win  = stft_params.win;
    target_freq = 1500; 
    [~, bin_idx] = min(abs(linspace(0, fs/2, nfft/2+1) - target_freq));

    % --- 1. Get Data Vector ---
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
    
    % --- 3. Test Reversed Shifts ---
    fprintf('Testing Reversed Permutations...\n');
    
    % The Base Reversed Order (Mic 1 stays at 1)
    % [1 6 5 4 3 2] is standard CW if Mic 1 is front.
    % But maybe Mic 1 is not front in the simulation? 
    % We test pure reversed vector: [6 5 4 3 2 1] and shifts.
    
    %base_rev = [7,6, 5, 4, 3, 2, 1];
    base_rev = [1,2,3,4,5,6,7];
    
    best_score = 0;
    best_perm = [];
    
    for s = 0:Q-1
        % Create permutation: Circular shift of the reversed base
        perm = circshift(base_rev, s);
        
        % Measure Match
        v_test = v_model(perm);
        score = abs(v_test' * v_data);
        
        fprintf('  Shift %d -> Indices %s -> Score: %.4f\n', s, mat2str(perm), score);
        
        if score > best_score
            best_score = score;
            best_perm = perm;
        end
    end
    
    fprintf('--------------------------------------------------\n');
    if best_score > 0.9
        fprintf('SOLVED! Use these indices: %s\n', mat2str(best_perm));
        calibrated_indices = best_perm;
    else
        fprintf('FAILURE: Best score is only %.3f. Check 30-degree rotation?\n', best_score);
        calibrated_indices = 1:Q;
    end
end