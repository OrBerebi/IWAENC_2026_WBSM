function [dist, az, el] = calculate_doa(array_pos, source_pos)
    % array_pos:  [x, y, z] of the array
    % source_pos: [x, y, z] of the source
    
    % 1. Calculate the relative vector (Source relative to Array)
    v = source_pos - array_pos;
    
    % 2. Calculate Distance (Euclidean Norm)
    dist = norm(v);
    
    % Check for singularity if source is at the same position as array
    if dist < 1e-9
        az = 0; el = 0;
        return;
    end
    
    % 3. Calculate Azimuth (phi)
    % Range: [0, 2*pi], 0 is Front (+x), CCW increase
    az = atan2(v(2), v(1));
    az = mod(az, 2*pi); 
    
    % 4. Calculate Elevation (theta)
    % Range: [0, pi], 0 is Up (+z), pi/2 is Horizontal, pi is Down (-z)
    el = acos(v(3) / dist);
end