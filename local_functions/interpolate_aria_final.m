function aria_lebedev_tmp = interpolate_aria_final(aria_atf)
% INTERPOLATE_ARIA_FINAL Robust SH Interpolation with Unit Correction.
%   Corrects Degree/Radian mismatches to fix the "Rubbish" SDR.
%   Target SDR: > 20 dB.

    %% 1. Configuration
    % Order 15 is extremely safe for 1674 points. 
    % It captures everything up to ~3-4 kHz perfectly and smooths the rest.
    SH_ORDER = 15; 
    REG_LAMBDA = 1e-3; % Standard regularization

    % Target Grid
    LEBEDEV_ORDER_OUT = 35; 

    fprintf('--- Starting Robust SH Interpolation (Order %d) ---\n', SH_ORDER);

    %% 2. Process Source Coordinates (With Unit Safety)
    src_az = aria_atf.sourceGrid.azimuth(:);   
    src_el = aria_atf.sourceGrid.elevation(:); % 0=Up, 180=Down
    src_r  = aria_atf.sourceGrid.r(:);

    % --- CRITICAL FIX: Detect and Fix Degrees ---
    if max(abs(src_az)) > 7.0 || max(abs(src_el)) > 4.0
        fprintf('-> Detected DEGREES in source grid. Converting to RADIANS...\n');
        src_az = deg2rad(src_az);
        src_el = deg2rad(src_el);
    else
        fprintf('-> Detected RADIANS in source grid.\n');
    end
    
    % Generate SH Basis for Source
    fprintf('Calculating SH Basis for Source Grid...\n');
    Y_src = get_SH_basis(SH_ORDER, src_az, src_el);
    
    % Compute Inverse (Encoder)
    % W_enc = (Y'Y + lambda*I)^-1 * Y'
    fprintf('Inverting Matrix (Reg: %.1e)...\n', REG_LAMBDA);
    Y_gram = Y_src' * Y_src;
    Y_pinv = (Y_gram + REG_LAMBDA * eye(size(Y_gram))) \ Y_src';

    %% 3. Generate Target Grid (Lebedev)
    fprintf('Generating Target Lebedev Grid (Order %d)...\n', LEBEDEV_ORDER_OUT);
    [gridData, ~, ~] = sofia_lebedev(Order_to_degree_lebedev(LEBEDEV_ORDER_OUT), 0);
    
    tgt_az = gridData(:,1); 
    tgt_el = gridData(:,2); % SOFiA usually returns Radians.
    
    % Generate SH Basis for Target
    Y_tgt = get_SH_basis(SH_ORDER, tgt_az, tgt_el);

    %% 4. Interpolate
    % Weights Matrix: [Target x Source]
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
    aria_lebedev.grid.elevation = tgt_el; 
    aria_lebedev.grid.r = ones(size(tgt_az));
    aria_lebedev.fs = aria_atf.fs;
    aria_lebedev.type = sprintf('SH_Order%d_Corrected', SH_ORDER);

    aria_lebedev_tmp = aria_atf;
    aria_lebedev_tmp.data = aria_lebedev.data;
    aria_lebedev_tmp.sourceGrid.azimuth = aria_lebedev.grid.azimuth;
    aria_lebedev_tmp.sourceGrid.elevation = aria_lebedev.grid.elevation;
    aria_lebedev_tmp.sourceGrid.r = aria_lebedev.grid.r;
    aria_lebedev_tmp.nData = length(aria_lebedev.grid.r);
end

%% --- HELPER: Real Spherical Harmonics ---
function Y = get_SH_basis(N, az, colat)
    % az, colat must be in RADIANS
    n_points = length(az);
    n_coeffs = (N+1)^2;
    Y = zeros(n_points, n_coeffs);
    
    idx = 1;
    for n = 0:N
        % Legendre Polynomials P_n(cos(theta))
        % Input to legendre must be cos(colat)
        P = legendre(n, cos(colat')); 
        
        for m = -n:n
            P_nm = P(abs(m)+1, :)';
            
            % Normalization
            norm = sqrt((2*n+1)/(4*pi) * factorial(n-abs(m))/factorial(n+abs(m)));
            
            % Real SH Basis Functions
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