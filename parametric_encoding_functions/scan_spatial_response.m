function scan_spatial_response(mic_signals, ATF_direct, stft_params, fs)
% SCAN_SPATIAL_RESPONSE
% Performs a Delay-and-Sum beamforming scan across 360 degrees.
% This visually reveals the geometry mismatch (Rotation vs Reflection).

    [nSamples, Q] = size(mic_signals);
    nfft = stft_params.nfft;
    win  = stft_params.win;
    nBins = nfft/2 + 1;
    
    % --- 1. Extract Signal Data Vector (at 1.5 kHz) ---
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
    v_data = E(:, end); % Principal component of the microphone signal
    v_data = v_data / norm(v_data);

    % --- 2. Scan ATF Library ---
    % We assume ATF_direct is [Q x nBins x nAngles]
    % or we need to find the angular dimension.
    [d1, d2, d3] = size(ATF_direct);
    
    % Heuristic: The dimension that isn't Q (6) or nBins (513) is Angles
    if d3 > Q && d3 > 10
        nAngles = d3;
        dim_angle = 3;
    elseif d1 > 10
        nAngles = d1;
        dim_angle = 1;
        % Permute to standard [Q x nBins x Angles]
        ATF_direct = permute(ATF_direct, [2, 3, 1]); % Just a guess, unlikely
    else
        error('Could not identify Angle dimension in ATF_direct (expecting >10 angles)');
    end
    
    scan_power = zeros(1, nAngles);
    
    fprintf('Scanning spatial response across %d ATF angles...\n', nAngles);
    
    for ang = 1:nAngles
        % Get Steering Vector for this angle
        v_steer = squeeze(ATF_direct(:, bin_idx, ang));
        v_steer = v_steer / norm(v_steer);
        
        % Beamforming Power: |w' * v_data|^2
        % (Matched Filter response)
        scan_power(ang) = abs(v_steer' * v_data)^2;
    end
    
    % --- 3. Plot ---
    % Normalize
    scan_power = scan_power / max(scan_power);
    
    % Assume angles are linearly spaced 0..360 for plotting
    angles_deg = linspace(0, 360, nAngles);
    
    figure('Color','w', 'Name', 'Spatial Scan Diagnosis');
    polarplot(deg2rad(angles_deg), scan_power, 'LineWidth', 2);
    title(sprintf('Array Response Pattern (Freq: %d Hz)', target_freq));
    rticklabels({'0 (Front?)', '30', '60', '90', '120', '150', '180', '210', '240', '270', '300', '330'});
    thetaticks(0:30:330);
    
    % Find Peak
    [~, max_idx] = max(scan_power);
    peak_angle = angles_deg(max_idx);
    
    fprintf('Peak Response detected at Index: %d (approx %.1f deg)\n', max_idx, peak_angle);
    fprintf('If actual source is at 0 deg, you have a rotation of %.1f deg.\n', peak_angle);
end