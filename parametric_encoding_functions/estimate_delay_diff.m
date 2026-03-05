function [test_brir_aligned, delay_diff] = estimate_delay_diff(ref_brir, test_brir, fs)
    % ESTIMATE_DELAY_DIFF Aligns test_brir to ref_brir based on Group Delay.
    % Returns the aligned signal and the calculated delay difference.

    gd_freq_lim = [0.4e3, 2.5e3];
    
    % 1. Calculate Group Delay for Test
    group_delays = group_delay_multi_channel(test_brir);
    f_vec_gd = linspace(0, fs/2, size(group_delays, 1)).';
    
    [~, f_gd_idx_low]  = min(abs(f_vec_gd - gd_freq_lim(1)));
    [~, f_gd_idx_high] = min(abs(f_vec_gd - gd_freq_lim(2)));
    
    TOA_gd_test_lr = mean(group_delays(f_gd_idx_low:f_gd_idx_high, :)); % in samples 
    
    % 2. Calculate Group Delay for Ref
    group_delays = group_delay_multi_channel(ref_brir);
    TOA_gd_ref_lr = mean(group_delays(f_gd_idx_low:f_gd_idx_high, :)); % in samples 
    
    % 3. Calculate Difference
    % If Ref is "later" (larger TOA), we need to delay Test (positive shift).
    delay_diff = TOA_gd_ref_lr - TOA_gd_test_lr; % 1x2 [left_diff, right_diff]
    
    % 4. Apply Sub-Sample Shift in Frequency Domain
    test_brir_aligned = zeros(size(test_brir));
    [N, n_ch] = size(test_brir);
    
    % Construct Frequency Vector k for Standard FFT (0 to 2pi)
    % Handles 0, positive freqs, and negative freqs correctly
    if mod(N, 2) == 0
        k = [0:N/2-1, -N/2:-1]';
    else
        k = [0:(N-1)/2, -(N-1)/2:-1]';
    end
    
    for ch = 1:n_ch
        shift = delay_diff(ch);
        
        % Transform to Frequency Domain
        X = fft(test_brir(:, ch));
        
        % Apply Linear Phase Shift: exp(-j * w * shift)
        % w = 2*pi*k/N
        phase_shift = exp(-1j * 2 * pi * k * shift / N);
        X_shifted = X .* phase_shift;
        
        % Transform back to Time Domain (real part to remove numerical noise)
        test_brir_aligned(:, ch) = real(ifft(X_shifted));
    end
end