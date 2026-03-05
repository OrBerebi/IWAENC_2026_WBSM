function W = buildWeightedW(omega_V, omega_sector, alpha, beta)
% buildWeightedW Constructs a diagonal weighting matrix for weighted LS
%
%   W = buildWeightedW(omega_V, omega_sector, alpha, beta)
%
% Inputs:
%   omega_V        - (N_dirs x 3) array of [elevation, azimuth, distance] for all directions
%   omega_sector   - (n_sector x 3) array of [elevation, azimuth, distance] for in-sector directions
%   alpha          - scalar weight for sector directions
%   beta           - scalar weight for out-of-sector directions
%
% Output:
%   W              - (N_dirs x N_dirs) sparse diagonal weight matrix
%
% Example:
%   W = buildWeightedW(omega_V, omega_sector, 1, 0.01);

% Validate inputs
if size(omega_V,2) ~= 3 || size(omega_sector,2) ~= 3
    error('omega_V and omega_sector must be N x 3 arrays');
end

omega_V      = round(omega_V(:, 1:2), 4);
omega_sector = round(omega_sector(:, 1:2), 4);

% Identify sector indices
[inSector, loc] = ismember(omega_V, omega_sector, 'rows');
if sum(inSector) ~= size(omega_sector,1)
    error('Number of matching sector directions (%d) does not equal expected (%d)',...
        sum(inSector), size(omega_sector,1));
end

% Build weight vector
N = size(omega_V,1);
weights = beta * ones(N,1);
weights(inSector) = alpha;

% Construct sparse diagonal weighting matrix
W = spdiags(weights, 0, N, N);
end
