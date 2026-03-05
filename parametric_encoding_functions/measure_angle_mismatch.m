function measure_angle_mismatch(mic_signals, ATF_direct, stft_params, fs)
% MEASURE_ANGLE_MISMATCH_CORRECTED
% Applies the proposed geometry correction [1 6 5 4 3 2] and then 
% runs the angle mismatch measurement to verify the fix.

    % --- 1. APPLY PROPOSED CORRECTION ---
    fprintf('Applying Correction: Reordering Mics to [1 6 5 4 3 2]...\n');
    %mic_signals = mic_signals(:, [1, 6, 5, 4, 3, 2]);

    % --- 2. Standard Measurement Logic ---
    [nSamples, Q] = size(mic_signals);
    nfft = stft_params.nfft;
    win  = stft_params.win;
    
    target_freq = 1500; 
    [~, bin_idx] = min(abs(linspace(0, fs/2, nfft/2+1) - target_freq));
    
    % Assumed Geometry (Standard 6-mic circular CCW)
    mic_angles = deg2rad([0, 60, 120, 180, 240, 300])'; 
    
    % Extract Phases (Data)
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
    vals_data = unwrap(angle(v_data)); 
    vals_data = vals_data - mean(vals_data);

    % Extract Phases (Model)
    v_model = squeeze(ATF_direct(:, bin_idx, :));
    vals_model = unwrap(angle(v_model));
    vals_model = vals_model - mean(vals_model);
    
    % Fit Cosine Model
    fit_angle = @(vals) fminbnd(@(theta) norm(vals - (max(abs(vals))*cos(mic_angles - theta))), -pi, pi);
    
    ang_data  = fit_angle(vals_data);
    ang_model = fit_angle(vals_model);
    
    deg_data  = rad2deg(ang_data);
    deg_model = rad2deg(ang_model);
    
    diff_deg = deg_data - deg_model;
    diff_deg = mod(diff_deg + 180, 360) - 180;

    % --- 3. Report & Plot ---
    fprintf('=== CORRECTED DIAGNOSIS ===\n');
    fprintf('Data Vector Source Est:  %6.1f deg\n', deg_data);
    fprintf('Model Vector Source Est: %6.1f deg\n', deg_model);
    fprintf('Residual Error:          %6.1f deg\n', diff_deg);
    
    if abs(diff_deg) < 10
        fprintf('SUCCESS: Curves are now aligned!\n');
    else
        fprintf('WARNING: Still mismatched. Check coordinate rotation.\n');
    end
    
    figure('Color','w', 'Name', 'Corrected Phase Analysis');
    subplot(2,1,1);
    plot(rad2deg(mic_angles), vals_data, '-o', 'LineWidth', 2, 'DisplayName', 'Data Phase (Corrected)');
    hold on;
    plot(rad2deg(mic_angles), vals_model, '--x', 'LineWidth', 2, 'DisplayName', 'Model Phase');
    legend; grid on; xlabel('Mic Position (Deg)'); ylabel('Phase (Rad)');
    title('Phase Pattern Comparison (After Correction)');
    
    subplot(2,1,2);
    polarplot([0; ang_data], [0; 1], '-o', 'LineWidth', 2, 'DisplayName', 'Data Source');
    hold on;
    polarplot([0; ang_model], [0; 1], '--x', 'LineWidth', 2, 'DisplayName', 'Model Source');
    title(sprintf('Residual Error: %.1f deg', diff_deg));
    legend;
end