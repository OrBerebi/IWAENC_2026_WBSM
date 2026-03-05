function [c_l_scaled, c_r_scaled, alpha_gain] = match_bsm_magls_gain(c_l, c_r, V_k, h_ref_data, freqs_sig, f_cut_magLS)
% match_bsm_magls_gain.m
% Calculates the RMS energy correction factor strictly for the MagLS 
% high-frequency region and smoothly crossfades the gain to prevent spectral shelves.

    % 1. Setup frequency bands exactly like your main script
    f_ratio = 2^0.25;
    f_low = f_cut_magLS / f_ratio;
    f_high = f_cut_magLS * f_ratio;
    
    [~, n_freqs] = size(c_l);
    n_dirs = size(h_ref_data, 1);
    
    % 2. Isolate purely MagLS high frequencies for the calculation
    magls_idx = freqs_sig > f_high;
    
    h_hat_l_high = zeros(n_dirs, sum(magls_idx));
    h_hat_r_high = zeros(n_dirs, sum(magls_idx));
    h_ref_l_high = h_ref_data(:, magls_idx, 1);
    h_ref_r_high = h_ref_data(:, magls_idx, 2);
    
    % Reconstruct only the high-frequency HRTFs
    freq_indices = find(magls_idx);
    for i = 1:length(freq_indices)
        f = freq_indices(i);
        V_k_curr = V_k(:, :, f);
        if ~any(isnan(V_k_curr), 'all')
            h_hat_l_high(:, i) = (c_l(:, f)' * V_k_curr).';
            h_hat_r_high(:, i) = (c_r(:, f)' * V_k_curr).';
        end
    end
    
    % 3. Calculate Phase-Agnostic Energy Ratio (MagLS Region Only)
    hat_flat = [h_hat_l_high(:); h_hat_r_high(:)];
    ref_flat = [h_ref_l_high(:); h_ref_r_high(:)];
    
    valid_idx = ~isnan(hat_flat) & ~isnan(ref_flat);
    hat_valid = hat_flat(valid_idx);
    ref_valid = ref_flat(valid_idx);
    
    energy_ref = sum(abs(ref_valid).^2);
    energy_hat = sum(abs(hat_valid).^2);
    
    alpha_gain = sqrt(energy_ref / (energy_hat + 1e-12));
    fprintf('Calculated MagLS-Region Gain Factor: %.4f\n', alpha_gain);
    
    % 4. Create a frequency-dependent gain curve
    %    1.0 for LS, alpha_gain for MagLS, smoothly interpolated in between
    gain_curve = ones(1, n_freqs);
    for f = 1:n_freqs
        current_freq = freqs_sig(f);
        if current_freq > f_high
            gain_curve(f) = alpha_gain;
        elseif current_freq > f_low
            % Crossfade the gain exactly mirroring your filter crossfade
            alpha_cross = (current_freq - f_low) / (f_high - f_low);
            gain_curve(f) = (1 - alpha_cross) * 1.0 + alpha_cross * alpha_gain;
        end
    end
    
    % 5. Apply the smooth gain curve across the coefficients
    % (MATLAB will implicitly expand the 1D gain_curve across the rows of c)
    c_l_scaled = c_l .* gain_curve;
    c_r_scaled = c_r .* gain_curve;
end