function filtered_omega = filterOmegaByGCD(omega, angle, GCD_radius, k)
%FILTEROMEGABYGCD Filters omega entries based on great circle distance (GCD)
%   Inputs:
%       omega       - [N x 3] matrix of [elevation, azimuth, radius] in degrees
%                     elevation: 0 = +Z, 90 = XY plane, 180 = -Z
%       angle       - [1 x 2] vector of [elevation, azimuth] in same convention
%       GCD_radius  - scalar threshold in degrees
%   Output:
%       filtered_omega - subset of omega within GCD_radius (in degrees) of the look direction

    % Convert user-style elevation to MATLAB-style elevation:
    % el_matlab = 90 - el_user
    omega_el_rad = pi/2 - omega(:,1);
    omega_az_rad = omega(:,2);
    r = omega(:,3);

    % Convert to Cartesian coordinates
    [x, y, z] = sph2cart(omega_az_rad, omega_el_rad, r);
    points_cartesian = [x, y, z];

    % Look direction conversion
    angle_el_rad = deg2rad(90 - angle(1));
    angle_az_rad = deg2rad(angle(2));
    [lx, ly, lz] = sph2cart(angle_az_rad, angle_el_rad, r(1));
    look_vec = [lx, ly, lz];

    % Build KDTree using Cartesian coordinates
    tree = KDTreeSearcher(points_cartesian);

    % Query k nearest neighbors (Euclidean distance, approximates small-angle GCD)
    [index, euclid_dist] = knnsearch(tree, look_vec, 'K', k);


    GCD = rad2deg(2 * asin(euclid_dist / (2 * r(1))));
    
    mask = GCD <= GCD_radius;
    
    filtered_omega = omega(index(mask),:);




    

    
    % % Normalize vectors (important if r ≠ 1)
    % points_cartesian = normalize(points_cartesian, 2);
    % look_vec = look_vec / norm(look_vec);
    % 
    % % Compute cosine of angle between each point and the look direction
    % cos_theta = points_cartesian * look_vec';
    % cos_theta = max(min(cos_theta, 1), -1); % clamp for safety
    % 
    % % Angular distance in radians
    % angular_distances = acos(cos_theta);
    % 
    % % Threshold in radians
    % GCD_radius_rad = deg2rad(GCD_radius);
    % 
    % % Apply filter
    % mask = angular_distances <= GCD_radius_rad;
    % filtered_omega = omega(mask, :);
end
