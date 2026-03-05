function idx = findClosestDirection(th_deg, ph_deg, Omega)
% findClosestDirection - Finds the index of the closest direction in Omega
%
% Syntax:
%   idx = findClosestDirection(th_deg, ph_deg, Omega)
%
% Inputs:
%   th_deg - Elevation in degrees (0 = top of sphere, 90 = horizon)
%   ph_deg - Azimuth in degrees (0 = x-axis, counterclockwise)
%   Omega  - [N x 2] matrix of [elevation, azimuth] in **radians**
%
% Output:
%   idx    - Index in Omega of the closest direction

    % Convert input angles to radians
    th = deg2rad(th_deg);
    ph = deg2rad(ph_deg);
    
    % Convert all directions to Cartesian coordinates on unit sphere
    % Input direction
    x0 = sin(th) * cos(ph);
    y0 = sin(th) * sin(ph);
    z0 = cos(th);
    dir0 = [x0, y0, z0];

    % Omega directions
    el = Omega(:, 1);
    az = Omega(:, 2);
    x = sin(el) .* cos(az);
    y = sin(el) .* sin(az);
    z = cos(el);
    dirs = [x, y, z];

    % Compute Euclidean distances
    dists = sqrt(sum((dirs - dir0).^2, 2));

    % Find index of minimum distance
    [~, idx] = min(dists);
end
