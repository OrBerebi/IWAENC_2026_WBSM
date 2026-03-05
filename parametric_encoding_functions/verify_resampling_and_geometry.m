function verify_resampling_and_geometry(mic_signals, ATF_direct, stft_params, fs)
% VERIFY_RESAMPLING_AND_GEOMETRY
% Diagnoses Resampling errors and Geometric alignment for the specific
% Semi-Circular Array provided.
%
% ARGS:
%   mic_signals: [N x 6] Time domain signal
%   ATF_direct:  [6 x 513] Frequency domain vector (NFFT=1024)
%   stft_params: Struct with .nfft=4096, .win, .hop
%   fs:          Sampling Rate (e.g. 48000)

    % --- GEOMETRY CORRECTION ---
    % The Simulation and ATF have opposite Azimuth signs.
    % We reverse the semi-circular array to align them.
    % reverse_idx = [6, 5, 4, 3, 2, 1];
    % ATF_direct = ATF_direct(reverse_idx, :); 
    % 
    % --- 1. Constants & Geometry ---
    % Array Azimuths (deg) provided by user
    mic_az = [295.7143, 321.4286, 347.1429, 12.8571, 38.5714, 64.2857];

    %mic_az = [334.4339, 3.9972, 24.2876, 34.6341, 330.3173, 267.1545,93.5441];
    
    % Unwrap negative angles for plotting continuity (-180 to 180)
    mic_az_plot = mic_az;
    mic_az_plot(mic_az_plot > 180) = mic_az_plot(mic_az_plot > 180) - 360;
    
    [nSamples, Q] = size(mic_signals);
    nfft_sys = stft_params.nfft; % 4096
    nBins_sys = nfft_sys/2 + 1;  % 2049
    
    % --- 2. Extract Data Phase (Ground Truth) ---
    % We average the covariance over the first 50 frames (Direct Sound)
    target_freq = 1500; % Check at 1.5kHz
    bin_sys = round(target_freq / (fs/nfft_sys)) + 1;
    
    Cx = zeros(Q, Q);
    mic_pad = [mic_signals; zeros(nfft_sys, Q)];
    nFrames = min(50, floor((nSamples-nfft_sys)/stft_params.hop)); 
    
    for i = 0:nFrames-1
        idx = i*stft_params.hop + 1;
        x_win = mic_pad(idx:idx+nfft_sys-1, :) .* stft_params.win;
        X_f = fft(x_win, nfft_sys, 1);
        x_k = X_f(bin_sys, 1:Q).'; 
        Cx = Cx + (x_k * x_k');
    end
    [E, ~] = eig(Cx);
    v_data = E(:, end); % Dominant eigenvector
    
    % Reference to Mic 4 (Closest to Source 11.4 deg)
    % We expect Mic 4 to be 0 phase (or leading).
    ref_mic_idx = 4; 
    phase_data = unwrap(angle(v_data ./ v_data(ref_mic_idx)));

    % --- 3. Process Model (ATF) with Resampling ---
    % ATF is 6x513. We need 6x2049.
    [Q_atf, nBins_atf] = size(ATF_direct);
    
    % METHOD A: DFT Resampling (Time Padding) - What 'perform_resampling' does
    % This assumes the 513 bins cover 0..Nyquist of the SAME Fs as the 2049 bins.
    
    % Replicate 'perform_resampling' logic locally for transparency
    N_in  = (nBins_atf - 1) * 2; % 1024
    N_out = (nBins_sys - 1) * 2; % 4096
    
    ATF_perm = permute(ATF_direct, [2, 1]); % [513 x 6]
    ATF_full = [ATF_perm; conj(ATF_perm(end-1:-1:2, :))]; % [1024 x 6]
    h_temp = real(ifft(ATF_full, N_in, 1));
    h_pad  = [h_temp; zeros(N_out - N_in, Q)]; % Pad to 4096
    ATF_new_full = fft(h_pad, N_out, 1);
    v_model_resampled = ATF_new_full(bin_sys, :).'; % Extract bin at 1.5kHz
    
    phase_model = unwrap(angle(v_model_resampled ./ v_model_resampled(ref_mic_idx)));

    % --- 4. Plot ---
    figure('Color','w', 'Name', 'Semi-Circular Array Diagnostics');
    
    subplot(2,1,1);
    plot(mic_az_plot, phase_data, '-o', 'LineWidth', 2, 'DisplayName', 'Data (Simulation)');
    hold on;
    plot(mic_az_plot, phase_model, '--x', 'LineWidth', 2, 'DisplayName', 'Model (Resampled)');
    xline(11.4, 'g--', 'Source DOA (11.4^\circ)');
    grid on;
    xlabel('Mic Azimuth (Deg)'); ylabel('Relative Phase (Rad)');
    title(sprintf('Phase Alignment at %.1f Hz (Ref: Mic 4)', target_freq));
    legend('Location','best');
    
    % Annotation
    dim = [.15 .6 .3 .3];
    str = {'Checklist:', ...
           '1. Data Curve Minimum should be near 11.4 deg (Mic 4).', ...
           '2. Model Curve should overlap Data.', ...
           '3. If Model is flatter -> Resampling is stretching frequency.'};
    %annotation('textbox',dim,'String',str,'FitBoxToText','on', 'BackgroundColor', 'w');

    % --- 5. Report Delays ---
    fprintf('=== DIAGNOSTICS ===\n');
    fprintf('Source: 11.4 deg. Closest Mic: 4 (12.9 deg).\n');
    fprintf('Data Phase @ Mic 4: %.3f (Should be near 0)\n', phase_data(4));
    fprintf('Model Phase @ Mic 4: %.3f\n', phase_model(4));
    
    % Check Slope (Group Delay proxy) between Mic 3 (-13 deg) and Mic 4 (13 deg)
    slope_data = phase_data(4) - phase_data(3);
    slope_model = phase_model(4) - phase_model(3);
    
    fprintf('Phase Diff (Mic 4 - Mic 3):\n');
    fprintf('  Data:  %.3f rad\n', slope_data);
    fprintf('  Model: %.3f rad\n', slope_model);
    fprintf('  Ratio (Model/Data): %.2f\n', slope_model/slope_data);
    
    if abs(slope_model/slope_data - 1) > 0.1
        fprintf('\nFAIL: Significant Mismatch!\n');
        if abs(slope_model/slope_data - 3) < 0.5
            fprintf('--> Ratio is ~3.0. Is ATF sampled at 16kHz but System at 48kHz?\n');
            fprintf('    (1024 bins @ 16k = 4096 bins @ 48k? No, check ratios.)\n');
        else
            fprintf('--> Likely Geometry or Resampling Error.\n');
        end
    else
        fprintf('\nPASS: Slopes match. Geometry seems correct.\n');
    end
end