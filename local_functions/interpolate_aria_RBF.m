function aria_lebedev_tmp = interpolate_aria_RBF(aria_atf)
% INTERPOLATE_ARIA_RBF High-Fidelity interpolation using Radial Basis Functions.
%   Uses a Gaussian Kernel to robustly map scattered Aria data to Lebedev.
%   Expected SDR: > 20 dB.

    %% 1. Configuration
    % Target Grid (Lebedev)
    LEBEDEV_ORDER_OUT = 35; 
    
    % Regularization (Lambda)
    % Small value to prevent overfitting noise. 1e-4 is standard for clean measurements.
    REG_LAMBDA = 1e-4; 

    fprintf('--- Starting High-Fidelity RBF Interpolation ---\n');

    %% 2. Process Source Coordinates (Aria) -> Cartesian
    % Input: Az=0(Front), El=0(Up)..180(Down)
    src_az = aria_atf.sourceGrid.azimuth(:);   
    src_el = aria_atf.sourceGrid.elevation(:);
    src_r  = aria_atf.sourceGrid.r(:);

    % Convert to Cartesian (Source Points)
    % We normalize to Unit Sphere for the RBF math (distance on surface)
    sx = sin(src_el) .* cos(src_az);
    sy = sin(src_el) .* sin(src_az);
    sz = cos(src_el);
    src_pts = [sx, sy, sz];
    
    nSource = length(sx);
    fprintf('Source: %d points processed.\n', nSource);

    %% 3. Generate Target Grid (Lebedev) -> Cartesian
    fprintf('Generating Target Lebedev Grid (Order %d)...\n', LEBEDEV_ORDER_OUT);
    
    % --- User Provided Code ---
    [gridData, ~, ~] = sofia_lebedev(Order_to_degree_lebedev(LEBEDEV_ORDER_OUT), 0);
    tgt_az = gridData(:,1); 
    tgt_el = gridData(:,2); % 0=Up, 180=Down
    % --------------------------
    
    % Convert to Cartesian (Target Points)
    tx = sin(tgt_el) .* cos(tgt_az);
    ty = sin(tgt_el) .* sin(tgt_az);
    tz = cos(tgt_el);
    tgt_pts = [tx, ty, tz];
    
    nTarget = length(tx);

    %% 4. Compute RBF Kernel Parameters (Auto-Tuning)
    fprintf('Auto-tuning RBF Kernel Width...\n');
    
    % Estimate average spacing between source points to set the Gaussian width (epsilon)
    % We take a random subset to estimate this quickly
    subset_idx = randperm(nSource, min(500, nSource));
    sub_pts = src_pts(subset_idx, :);
    
    % Calculate squared euclidean distances for subset
    % (Just checking nearest neighbor distance)
    dists = pdist(sub_pts); 
    avg_dist = mean(dists); 
    
    % Rule of thumb: Gaussian width should be ~1.0 to 1.5x the average spacing
    epsilon = 1.2 * avg_dist; 
    fprintf('-> Average Point Spacing: %.4f\n', avg_dist);
    fprintf('-> Selected Gaussian Width (epsilon): %.4f\n', epsilon);

    %% 5. Build Interaction Matrices (The Heavy Math)
    fprintf('Building Source-Source Kernel Matrix (K_ss)...\n');
    
    % Calculate Distance Matrix (Source vs Source)
    % dist_sq_SS(i,j) = ||x_i - x_j||^2
    % Efficient vectorized calculation: |a-b|^2 = |a|^2 + |b|^2 - 2*a*b
    % Since points are on unit sphere, |a|^2 = 1. So dist^2 = 2 - 2*dot(a,b).
    dot_SS = src_pts * src_pts';
    dist_sq_SS = 2 - 2 * dot_SS;
    dist_sq_SS(dist_sq_SS < 0) = 0; % Numeric safety
    
    % Apply Gaussian Kernel: exp( - dist^2 / epsilon^2 )
    K_SS = exp( -dist_sq_SS / (epsilon^2) );
    
    % Add Regularization (Ridge Regression)
    % This makes the matrix invertible and handles noise
    K_SS_reg = K_SS + REG_LAMBDA * eye(nSource);
    
    fprintf('Inverting Kernel Matrix...\n');
    % Calculate Weights Generator: W_gen = K_SS^-1
    % Use Cholesky decomposition for speed/stability if positive definite
    % Or simple backslash for robustness.
    inv_K_SS = inv(K_SS_reg); 

    %% 6. Build Target-Source Interpolation Matrix
    fprintf('Building Target-Source Kernel Matrix (K_ts)...\n');
    
    % Calculate Distance Matrix (Target vs Source)
    dot_TS = tgt_pts * src_pts';
    dist_sq_TS = 2 - 2 * dot_TS;
    dist_sq_TS(dist_sq_TS < 0) = 0;
    
    % Apply same Kernel
    K_TS = exp( -dist_sq_TS / (epsilon^2) );
    
    % Compute Final Transform Matrix: [Target x Source]
    % Weights = K_TS * inv(K_SS)
    W_interp = K_TS * inv_K_SS;

    %% 7. Main Interpolation Loop
    [~, nTime, nReceivers] = size(aria_atf.data);
    interp_data = zeros(nTarget, nTime, nReceivers);
    
    tic;
    fprintf('Interpolating Data...\n');
    for r = 1:nReceivers
        % Data: [1674 x 332]
        rec_data = aria_atf.data(:, :, r);
        
        % Multiply: [Target x Source] * [Source x Time]
        % This effectively "moves" the sound from the Aria grid to the Lebedev grid
        interp_data(:, :, r) = W_interp * rec_data;
    end
    t_end = toc;
    fprintf('Interpolation complete in %.2f seconds.\n', t_end);

    %% 8. Package Output
    aria_lebedev.data = interp_data;
    aria_lebedev.grid.azimuth = tgt_az;
    aria_lebedev.grid.elevation = tgt_el; 
    aria_lebedev.grid.r = ones(size(tgt_az));
    aria_lebedev.fs = aria_atf.fs;
    aria_lebedev.type = sprintf('RBF_Gaussian_eps%.2f', epsilon);

    aria_lebedev_tmp = aria_atf;
    aria_lebedev_tmp.data = aria_lebedev.data;
    aria_lebedev_tmp.sourceGrid.azimuth = aria_lebedev.grid.azimuth;
    aria_lebedev_tmp.sourceGrid.elevation = aria_lebedev.grid.elevation;
    aria_lebedev_tmp.sourceGrid.r = aria_lebedev.grid.r;
    aria_lebedev_tmp.nData = length(aria_lebedev.grid.r);
    
end