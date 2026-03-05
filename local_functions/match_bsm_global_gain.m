function [c_l_scaled, c_r_scaled, alpha_gain] = match_bsm_global_gain(c_l, c_r, V_k, h_ref_data)
% match_bsm_global_gain.m
% Reconstructs the estimated HRTFs across all frequencies, calculates the 
% single optimal gain factor based on TOTAL ENERGY (RMS), and applies it.
% Safe for MagLS because it ignores phase mismatches.

    [~, n_freqs] = size(c_l);
    n_dirs = size(h_ref_data, 1);
    
    h_hat_l = zeros(n_dirs, n_freqs);
    h_hat_r = zeros(n_dirs, n_freqs);
    
    % 1. Reconstruct the estimated HRTFs
    for f = 1:n_freqs
        V_k_curr = V_k(:, :, f);
        if ~any(isnan(V_k_curr), 'all')
            h_hat_l(:, f) = (c_l(:, f)' * V_k_curr).';
            h_hat_r(:, f) = (c_r(:, f)' * V_k_curr).';
        end
    end
    
    % 2. Flatten both the estimate and the reference
    hat_flat = [h_hat_l(:); h_hat_r(:)];
    
    h_ref_l = h_ref_data(:, :, 1);
    h_ref_r = h_ref_data(:, :, 2);
    ref_flat = [h_ref_l(:); h_ref_r(:)];
    
    % 3. Filter out NaNs
    valid_idx = ~isnan(hat_flat) & ~isnan(ref_flat);
    hat_valid = hat_flat(valid_idx);
    ref_valid = ref_flat(valid_idx);
    
    % 4. Calculate Phase-Agnostic Energy Ratio
    % We match the total power of the estimate to the reference
    energy_ref = sum(abs(ref_valid).^2);
    energy_hat = sum(abs(hat_valid).^2);
    
    % alpha is the square root of the power ratio (RMS scaling)
    alpha_gain = sqrt(energy_ref / (energy_hat + 1e-12));
    
    fprintf('Calculated Energy-Matched Global Gain: %.4f\n', alpha_gain);
    
    % 5. Bake the gain directly into the coefficients
    c_l_scaled = alpha_gain * c_l;
    c_r_scaled = alpha_gain * c_r;
end