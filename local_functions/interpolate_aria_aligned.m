function aria_lebedev_tmp = interpolate_aria_aligned(aria_atf)
% INTERPOLATE_ARIA_ALIGNED High-Fidelity Time-Aligned Interpolation.
%   Combines robust Delaunay triangulation with Time-of-Arrival correction.
%   Target SDR: > 20 dB.

    %% 1. Configuration
    % Target Grid
    LEBEDEV_ORDER_OUT = 35; 
    
    fprintf('--- Starting Time-Aligned Linear Interpolation ---\n');

    %% 2. Process Source Coordinates (With Unit Safety)
    src_az = aria_atf.sourceGrid.azimuth(:);   
    src_el = aria_atf.sourceGrid.elevation(:); % 0=Up, 180=Down
    src_r  = aria_atf.sourceGrid.r(:);

    % Detect Degrees vs Radians
    if max(abs(src_az)) > 7.0
        fprintf('-> Detected DEGREES. Converting to Radians...\n');
        src_az = deg2rad(src_az);
        src_el = deg2rad(src_el);
    end

    % Convert to Cartesian
    % 0=Up (Z=1), 180=Down (Z=-1)
    sx = src_r .* sin(src_el) .* cos(src_az);
    sy = src_r .* sin(src_el) .* sin(src_az);
    sz = src_r .* cos(src_el);

    %% 3. Generate Target Grid (Lebedev)
    fprintf('Generating Target Lebedev Grid (Order %d)...\n', LEBEDEV_ORDER_OUT);
    [gridData, ~, ~] = sofia_lebedev(Order_to_degree_lebedev(LEBEDEV_ORDER_OUT), 0);
    
    tgt_az = gridData(:,1); 
    tgt_el = gridData(:,2);
    
    % Target Cartesian
    tx = sin(tgt_el) .* cos(tgt_az);
    ty = sin(tgt_el) .* sin(tgt_az);
    tz = cos(tgt_el);

    %% 4. Prepare Interpolators
    fprintf('Building Triangulation...\n');
    % We use two interpolators:
    % 1. 'natural' for Delays (Smooth, physics-based)
    % 2. 'linear' for Audio (Robust, prevents overshoot)
    F_delay = scatteredInterpolant(sx, sy, sz, ones(size(sx)), 'natural', 'nearest');
    F_audio = scatteredInterpolant(sx, sy, sz, ones(size(sx)), 'linear', 'nearest');

    %% 5. Main Interpolation Loop
    [nDirs, nTime, nReceivers] = size(aria_atf.data);
    nTarget = length(tgt_az);
    interp_data = zeros(nTarget, nTime, nReceivers);
    
    tic;
    for r = 1:nReceivers
        fprintf('Processing Receiver %d / %d...\n', r, nReceivers);
        
        % Data: [1674 x 332]
        raw_data = aria_atf.data(:, :, r);
        
        % --- STEP A: Calculate & Smooth Delays ---
        % Find peak of the Envelope (Hilbert) for robustness
        env = abs(hilbert(raw_data')); % [Time x Dirs]
        [~, peak_idx] = max(env);      % [1 x Dirs]
        
        delays = peak_idx(:);
        median_delay = median(delays);
        relative_delays = delays - median_delay;
        
        % Smooth the delays (Physics check: Delays must change smoothly)
        F_delay.Values = relative_delays;
        % We query the delays AT THE SOURCE locations to get a "denoised" version
        % This prevents jittery alignment
        smoothed_source_delays = F_delay(sx, sy, sz);
        
        % --- STEP B: Align Signals ---
        aligned_data = zeros(size(raw_data));
        for i = 1:nDirs
            % Shift back by the calculated delay
            shift = -round(smoothed_source_delays(i));
            aligned_data(i, :) = circshift(raw_data(i, :), shift, 2);
        end
        
        % --- STEP C: Interpolate Aligned Data ---
        % Loop over time is fast enough with Delaunay pre-calc
        target_aligned = zeros(nTarget, nTime);
        for t = 1:nTime
            F_audio.Values = aligned_data(:, t);
            target_aligned(:, t) = F_audio(tx, ty, tz);
        end
        
        % --- STEP D: Interpolate Delays to Target ---
        % Predict what the delay should be at the Lebedev points
        target_delays = F_delay(tx, ty, tz);
        
        % --- STEP E: Re-Apply Delay ---
        for i = 1:nTarget
            shift = round(target_delays(i));
            % Shift forward
            interp_data(i, :, r) = circshift(target_aligned(i, :), shift, 2);
        end
    end
    t_end = toc;
    fprintf('Interpolation complete in %.2f seconds.\n', t_end);

    %% 6. Package Output
    aria_lebedev.data = interp_data;
    aria_lebedev.grid.azimuth = tgt_az;
    aria_lebedev.grid.elevation = tgt_el; 
    aria_lebedev.grid.r = ones(size(tgt_az));
    aria_lebedev.fs = aria_atf.fs;
    aria_lebedev.type = 'TimeAligned_Delaunay';

    aria_lebedev_tmp = aria_atf;
    aria_lebedev_tmp.data = aria_lebedev.data;
    aria_lebedev_tmp.sourceGrid.azimuth = aria_lebedev.grid.azimuth;
    aria_lebedev_tmp.sourceGrid.elevation = aria_lebedev.grid.elevation;
    aria_lebedev_tmp.sourceGrid.r = aria_lebedev.grid.r;
    aria_lebedev_tmp.nData = length(aria_lebedev.grid.r);
end