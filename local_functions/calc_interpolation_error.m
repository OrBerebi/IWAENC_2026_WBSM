function [nmse_db, sdr] = calc_interpolation_error(input_struct)
% CALC_INTERPOLATION_ERROR Robust estimation of interpolation error.
%   Works with both raw 'aria_atf' and interpolated 'aria_lebedev' structs.

    fprintf('--- Starting Error Estimation (10%% Hold-Out) ---\n');

    % %% 1. Universal Field Detection (Fixes struct mismatch)
    % if isfield(input_struct, 'sourceGrid')
    %     % Case A: Original Aria Data
    %     grid_struct = input_struct.sourceGrid;
    %     fprintf('-> Detected Raw Data structure (.sourceGrid)\n');
    % elseif isfield(input_struct, 'grid')
    %     % Case B: Interpolated Lebedev Data
    %     grid_struct = input_struct.grid;
    %     fprintf('-> Detected Interpolated structure (.grid)\n');
    % else
    %     error('Input struct has neither .sourceGrid nor .grid fields.');
    % end

    %% 2. Extract and Sanitize Coordinates
    % We force everything to COLUMN vectors (:) to prevent size mismatch errors
    az = input_struct.sourceGrid.azimuth(:);
    el = input_struct.sourceGrid.elevation(:);
    
    % % handle Radius (scalar or vector)
    % if isfield(grid_struct, 'r')
    % 
    %     if length(r) == 1, r = repmat(r, length(az), 1); end
    % else
    %     r = ones(size(az)); % Default to unit sphere
    % end

    r = input_struct.sourceGrid.r(:);

    % Convert to Cartesian (0=Up/Colatitude logic)
    % Z = r * cos(theta)
    z = r .* cos(el);
    % R_xy = r * sin(theta)
    r_xy = r .* sin(el);
    % X, Y
    x = r_xy .* cos(az);
    y = r_xy .* sin(az);

    %% 3. Verify Data Dimensions
    [nDirs, nTime, nReceivers] = size(input_struct.data);
    
    if length(x) ~= nDirs
        error('Mismatch! Grid has %d points but Data has %d directions.', length(x), nDirs);
    end

    %% 4. Create Random Split
    holdout_percentage = 0.10;
    nTest = round(nDirs * holdout_percentage);
    
    rand_indices = randperm(nDirs);
    
    % Force indices to be COLUMNS
    idx_test  = rand_indices(1:nTest)';
    idx_train = rand_indices(nTest+1:end)';
    
    fprintf('Splitting Grid: %d Training points, %d Test points.\n', length(idx_train), length(idx_test));

    %% 5. Build Interpolator
    % We create the interpolator ONE time using the Training locations.
    % We pass dummy values (ones) just to initialize the structure.
    
    % CRITICAL FIX: Ensure the dummy values vector matches the x(idx_train) length/shape exactly.
    train_x = x(idx_train);
    train_y = y(idx_train);
    train_z = z(idx_train);
    dummy_vals = ones(size(train_x)); 

    % 'linear' is standard, 'nearest' prevents NaN at edges
    F = scatteredInterpolant(train_x, train_y, train_z, dummy_vals, 'linear', 'nearest');

    %% 6. Validation Loop
    total_error_energy = 0;
    total_true_energy  = 0;
    
    % Pre-fetch test coordinates
    test_x = x(idx_test);
    test_y = y(idx_test);
    test_z = z(idx_test);
    
    for r = 1:nReceivers
        rec_data = input_struct.data(:, :, r);
        
        for t = 1:nTime
            % 1. Train: Update the interpolator with the KNOWN values (Training set)
            % Force column vector to match the coordinates
            vals_train = rec_data(idx_train, t); 
            F.Values = vals_train;
            
            % 2. Predict: Estimate values at the HIDDEN locations (Test set)
            vals_predicted = F(test_x, test_y, test_z);
            
            % 3. Compare: Get actual hidden values
            vals_actual = rec_data(idx_test, t);
            
            % 4. Accumulate Error
            error_diff = vals_actual - vals_predicted;
            
            total_error_energy = total_error_energy + sum(abs(error_diff).^2);
            total_true_energy  = total_true_energy  + sum(abs(vals_actual).^2);
        end
    end

    %% 7. Results
    nmse_linear = total_error_energy / total_true_energy;
    nmse_db = 10 * log10(nmse_linear);
    sdr = -nmse_db;

    fprintf('------------------------------------------------\n');
    fprintf('Interpolation Error (Self-Consistency Test):\n');
    fprintf('NMSE: %.2f dB\n', nmse_db);
    fprintf('SDR:  %.2f dB\n', sdr);
    fprintf('------------------------------------------------\n');
end