function analyze_atf_mismatch(mic_signals, ATF_direct, stft_params, fs)
% ANALYZE_ATF_MISMATCH
% Compares the Theoretical Steering Vector (ATF_direct) against the 
% Actual Data-Derived Steering Vector (Principal Eigenvector).
%
% USAGE:
%   analyze_atf_mismatch(mic_signals, ATF_direct, stft_params, 16000)


    % --- GEOMETRY CORRECTION ---
    % The Simulation and ATF have opposite Azimuth signs.
    % We reverse the ATF to align them.
    %reverse_idx = [6, 5, 4, 3, 2, 1];
    %ATF_direct = ATF_direct(reverse_idx, :); 

    [nSamples, Q] = size(mic_signals);
    nfft = stft_params.nfft;
    hop  = stft_params.hop;
    win  = stft_params.win;
    nBins = nfft/2 + 1;
    freqs = linspace(0, fs/2, nBins);

    % --- 1. Resample ATF if needed ---
    [Q_v, nBins_v, ~] = size(ATF_direct);
    if nBins_v ~= nBins
        ATF_direct = perform_resampling(ATF_direct, nBins_v, nBins);
    end

    % --- 2. Estimate Data Covariance (Averaged) ---
    Cx_sum = zeros(Q, Q, nBins);
    
    % Only process first 0.5 seconds (Direct Sound dominant)
    limit_frames = min(50, floor((nSamples-nfft)/hop)); 
    
    mic_pad = [mic_signals; zeros(nfft, Q)];
    
    for i = 0:limit_frames-1
        idx = i*hop + 1;
        x_win = mic_pad(idx:idx+nfft-1, :) .* win;
        X_f = fft(x_win, nfft, 1);
        X_half = X_f(1:nBins, :).';
        
        for k = 1:nBins
            x_k = X_half(:, k);
            Cx_sum(:, :, k) = Cx_sum(:, :, k) + (x_k * x_k');
        end
    end
    
    % --- 3. Compare Vectors ---
    alignment_score = zeros(1, nBins);
    
    fprintf('Comparing ATF_direct vs Data...\n');
    
    for k = 2:nBins % Skip DC
        R = Cx_sum(:, :, k);
        v_model = squeeze(ATF_direct(:, k, :));
        
        % Data Vector (Principal Eigenvector = Signal Direction)
        [E, D] = eig(R);
        [~, max_idx] = max(diag(D));
        v_data = E(:, max_idx);
        
        % Normalize
        v_model = v_model / norm(v_model);
        v_data  = v_data  / norm(v_data);
        
        % Alignment Score (Dot Product)
        % 1.0 = Perfect Match
        % 0.0 = Orthogonal (Complete Mismatch)
        score = abs(v_model' * v_data);
        alignment_score(k) = score(1);
    end
    
    % --- 4. Plot ---
    figure('Name', 'Steering Vector Diagnostic', 'Color', 'w');
    plot(freqs, alignment_score, 'LineWidth', 2);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Alignment Score (0 to 1)');
    title('ATF Mismatch Analysis');
    subtitle('1.0 = Perfect, < 0.8 = Audible Degradation, < 0.2 = Signal Nulling');
    ylim([0 1.1]);
    yline(0.9, 'g--', 'Excellent');
    yline(0.5, 'r--', 'Severe Mismatch');
end

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