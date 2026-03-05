function aria_lebedev_tmp = interpolate_aria(aria_atf)
% INTERPOLATE_ARIA Interpolates Aria ATF data to a Lebedev grid using SOFiA.
%
%   INPUT:
%   aria_atf struct with:
%       .data : [1674 x 332 x 7] (Directions x Time x Receivers)
%       .sourceGrid.azimuth   : [1674 x 1] Radians (0=Front, 90=Left)
%       .sourceGrid.elevation : [1674 x 1] Radians (0=Up, 90=Front, 180=Down)
%       .sourceGrid.r         : [1674 x 1] Meters
%       .fs                   : Scalar (e.g. 48000)
%
%   OUTPUT:
%   aria_lebedev struct

    fprintf('--- Starting Interpolation (Lebedev Order 35) ---\n');

    %% 1. Process Source Coordinates (Aria Input)
    % Input format: 
    % Azimuth: 0=Front, 90=Left (Radians)
    % Elevation: 0=Up, 90=Front, 180=Down (Colatitude/Theta)
    
    src_az = aria_atf.sourceGrid.azimuth;   
    src_el = aria_atf.sourceGrid.elevation; % Colatitude
    src_r  = aria_atf.sourceGrid.r;

    % Convert Source Spherical (Colatitude) -> Cartesian
    % z = r * cos(theta)
    src.z = src_r .* cos(src_el);
    
    % x = r * sin(theta) * cos(phi)
    src.x = src_r .* sin(src_el) .* cos(src_az);
    
    % y = r * sin(theta) * sin(phi)
    src.y = src_r .* sin(src_el) .* sin(src_az);
    
    fprintf('Source Grid: %d points processed.\n', length(src.x));

    %% 2. Generate Target Grid (Lebedev via SOFiA)
    % --- User Provided Code Start ---
    lebedev_order = 35;
    % Note: Ensure sofia_lebedev and Order_to_degree_lebedev are in your path
    [gridData_Inf, ~, ~] = sofia_lebedev(Order_to_degree_lebedev(lebedev_order), 0);
    
    ph_lebedev = gridData_Inf(:,1); % Azimuth (Radians)
    th_lebedev = gridData_Inf(:,2); % Colatitude/Theta (Radians)
    % --- User Provided Code End ---

    fprintf('Target Grid: Generated %d Lebedev points.\n', length(ph_lebedev));

    % Convert Target Spherical (Colatitude) -> Cartesian (Unit Sphere)
    tgt.z = cos(th_lebedev);
    tgt.x = sin(th_lebedev) .* cos(ph_lebedev);
    tgt.y = sin(th_lebedev) .* sin(ph_lebedev);

    %% 3. Prepare Interpolator
    
    % 'linear' interpolation is standard.
    % 'nearest' extrapolation handles any tiny precision errors at the convex hull boundaries.
    F = scatteredInterpolant(src.x, src.y, src.z, ones(size(src.x)), 'linear', 'nearest');

    %% 4. Main Interpolation Loop
    [nDirs, nTime, nReceivers] = size(aria_atf.data);
    nTarget = length(tgt.x);
    
    % Pre-allocate Output: [TargetPoints x Time x Receivers]
    interp_data = zeros(nTarget, nTime, nReceivers);
    
    tic;
    for r = 1:nReceivers
        % Extract data for this receiver: [1674 x 332]
        rec_data = aria_atf.data(:, :, r);
        
        % Loop through time samples
        for t = 1:nTime
            % Update interpolator values
            F.Values = rec_data(:, t);
            
            % Interpolate onto Lebedev grid
            interp_data(:, t, r) = F(tgt.x, tgt.y, tgt.z);
        end
        
        fprintf('Receiver %d / %d processed.\n', r, nReceivers);
    end
    t_end = toc;
    fprintf('Interpolation complete in %.2f seconds.\n', t_end);

    %% 5. Package Output
    aria_lebedev.data = interp_data;
    
    % Save grid in the same format as input (Azimuth + Colatitude)
    aria_lebedev.grid.azimuth = ph_lebedev;
    aria_lebedev.grid.elevation = th_lebedev; % 0=Up, 180=Down
    aria_lebedev.grid.r = ones(size(ph_lebedev));
    
    aria_lebedev.fs = aria_atf.fs;
    aria_lebedev.type = sprintf('Lebedev_Order_%d', lebedev_order);

    aria_lebedev_tmp = aria_atf;
    aria_lebedev_tmp.data = aria_lebedev.data;
    aria_lebedev_tmp.sourceGrid.azimuth = aria_lebedev.grid.azimuth;
    aria_lebedev_tmp.sourceGrid.elevation = aria_lebedev.grid.elevation;
    aria_lebedev_tmp.sourceGrid.r = aria_lebedev.grid.r;
    aria_lebedev_tmp.nData = length(aria_lebedev.grid.r);

end