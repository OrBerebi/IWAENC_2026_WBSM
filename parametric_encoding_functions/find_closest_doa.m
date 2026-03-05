function idx = find_closest_doa(matrix_doa, x)
    % matrix_doa: [N x 2] matrix where [elevation, azimuth] are in rows
    % x: [1 x 2] row vector [elevation, azimuth]
    % Returns: Index of the row in matrix_doa closest to x
    
    % Extract components for clarity
    el_grid = matrix_doa(:, 1);
    az_grid = matrix_doa(:, 2);
    
    el_target = x(1);
    az_target = x(2);
    
    % 1. Convert Grid and Target to Unit Vectors (Cartesian)
    % Using your convention: el=0 is +z, az=0 is +x
    [grid_x, grid_y, grid_z] = sph_to_cart_custom(el_grid, az_grid);
    [target_x, target_y, target_z] = sph_to_cart_custom(el_target, az_target);
    
    % 2. Calculate Dot Product (Projection)
    % Since they are unit vectors, dot product = cos(alpha)
    % A larger dot product means a smaller angle/distance
    dot_products = (grid_x * target_x) + (grid_y * target_y) + (grid_z * target_z);
    
    % 3. Find the maximum similarity (minimum angular distance)
    [~, idx] = max(dot_products);
end

function [vx, vy, vz] = sph_to_cart_custom(el, az)
    % Helper to convert your specific convention to Cartesian unit vectors
    % Elevation (theta): 0 is Up (+z)
    % Azimuth (phi): 0 is Front (+x)
    vz = cos(el);
    vx = sin(el) .* cos(az);
    vy = sin(el) .* sin(az);
end