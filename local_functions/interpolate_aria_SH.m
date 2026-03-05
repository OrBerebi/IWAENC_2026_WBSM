function aria_lebedev_tmp = interpolate_aria_SH(aria_atf)
% INTERPOLATE_ARIA_SH_STABLE Robust SH interpolation.
%   Fixes the instability by reducing SH Order and adding regularization.

    %% 1. Configuration (The Safety Settings)
    % Order 35 was suicide. Order 20 is the "Goldilocks" zone for 1600 points.
    SH_ORDER_FIT = 20;  
    
    % Regularization (Lambda). 
    % 1e-6 is too weak. 1e-2 is standard for measured data.
    REG_LAMBDA = 1e-2;  

    % Target Grid (Lebedev)
    LEBEDEV_ORDER_OUT = 35; 

    fprintf('--- Starting Robust SH Interpolation (Order %d) ---\n', SH_ORDER_FIT);

    %% 2. Process Source Coordinates
    src_az = aria_atf.sourceGrid.azimuth(:);   
    src_colat = aria_atf.sourceGrid.elevation(:); % 0=Up, 180=Down
    
    % Generate SH Basis for Source
    % Matrix Size: [1674 x 441]
    fprintf('Calculating SH Basis for Source Grid...\n');
    Y_src = get_SH_basis(SH_ORDER_FIT, src_az, src_colat);
    
    % Calculate Pseudo-Inverse with Stronger Regularization
    % W_encode = (Y'Y + lambda*I)^-1 * Y'
    fprintf('Inverting Matrix (Regularization: %.1e)...\n', REG_LAMBDA);
    Y_gram = Y_src' * Y_src;
    Y_pinv = (Y_gram + REG_LAMBDA * eye(size(Y_gram))) \ Y_src';

    %% 3. Generate Target Grid (Lebedev)
    fprintf('Generating Target Lebedev Grid (Order %d)...\n', LEBEDEV_ORDER_OUT);
    % User provided code
    [gridData, ~, ~] = sofia_lebedev(Order_to_degree_lebedev(LEBEDEV_ORDER_OUT), 0);
    tgt_az = gridData(:,1); 
    tgt_colat = gridData(:,2);
    
    % Generate SH Basis for Target
    Y_tgt = get_SH_basis(SH_ORDER_FIT, tgt_az, tgt_colat);

    %% 4. Main Interpolation
    % Transform Matrix: [Target x Source]
    % This maps raw microphone signals directly to Lebedev pixels
    W_interp = Y_tgt * Y_pinv; 
    
    [nDirs, nTime, nReceivers] = size(aria_atf.data);
    nTarget = length(tgt_az);
    interp_data = zeros(nTarget, nTime, nReceivers);
    
    tic;
    for r = 1:nReceivers
        % Data: [1674 x 332]
        rec_data = aria_atf.data(:, :, r);
        
        % Multiply: [Target x 1674] * [1674 x Time]
        interp_data(:, :, r) = W_interp * rec_data;
    end
    t_end = toc;
    fprintf('Interpolation complete in %.2f seconds.\n', t_end);

    %% 5. Package Output
    aria_lebedev.data = interp_data;
    aria_lebedev.grid.azimuth = tgt_az;
    aria_lebedev.grid.elevation = tgt_colat; 
    aria_lebedev.grid.r = ones(size(tgt_az));
    aria_lebedev.fs = aria_atf.fs;
    aria_lebedev.type = sprintf('SH_Order%d_Reg%.2e', SH_ORDER_FIT, REG_LAMBDA);

    aria_lebedev_tmp = aria_atf;
    aria_lebedev_tmp.data = aria_lebedev.data;
    aria_lebedev_tmp.sourceGrid.azimuth = aria_lebedev.grid.azimuth;
    aria_lebedev_tmp.sourceGrid.elevation = aria_lebedev.grid.elevation;
    aria_lebedev_tmp.sourceGrid.r = aria_lebedev.grid.r;
    aria_lebedev_tmp.nData = length(aria_lebedev.grid.r);
end

%% --- HELPER: Real SH Basis ---
function Y = get_SH_basis(N, az, colat)
    n_points = length(az);
    n_coeffs = (N+1)^2;
    Y = zeros(n_points, n_coeffs);
    
    idx = 1;
    for n = 0:N
        P = legendre(n, cos(colat'));
        for m = -n:n
            P_nm = P(abs(m)+1, :)';
            norm = sqrt((2*n+1)/(4*pi) * factorial(n-abs(m))/factorial(n+abs(m)));
            
            if m > 0
                trig = cos(m * az);
                norm = norm * sqrt(2);
            elseif m < 0
                trig = sin(abs(m) * az);
                norm = norm * sqrt(2);
            else
                trig = ones(size(az));
            end
            Y(:, idx) = norm .* P_nm .* trig;
            idx = idx + 1;
        end
    end
end