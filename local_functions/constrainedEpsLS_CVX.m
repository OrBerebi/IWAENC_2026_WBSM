function C = constrainedEpsLS_CVX(V, h, omega_V, omega_sector, target_nmse_db)
% constrainedEpsLS_CVX Solves:
%   min_c ||V_s c - conj(h_s)||^2
%   s.t.  ||V_o c - conj(h_o)||^2 <= epsilon
%
% Inputs:
%   V             - [M x N x F] steering matrix
%   h             - [N x F] complex HRTF response
%   omega_V       - [N x 3] all directions
%   omega_sector  - [S x 3] directions inside the sector
%   epsilon       - scalar, energy bound for out-of-sector error
%
% Output:
%   C             - [M x F] filter coefficients

% Make sure CVX is available
assert(exist('cvx_begin', 'file') == 2, 'CVX not found. Add it to your path.');

% Find which directions are in/out of sector
[inSector, ~] = ismember(omega_V, omega_sector, 'rows');
sector_idx    = find(inSector);
outsector_idx = find(~inSector);

% Sizes
[M, N, F] = size(V);
C = zeros(M, F);  % output filters


N_o = length(outsector_idx);


% Solve for each frequency
for f = 1:F
    cvx_clear

    Vf = squeeze(V(:, :, f));   % M x N
    hf = conj(h(:, f));         % N x 1

    % Split matrices
    Vs = Vf(:, sector_idx);     % M x S
    hs = hf(sector_idx);        % S x 1

    Vo = Vf(:, outsector_idx);  % M x O
    ho = hf(outsector_idx);     % O x 1

    %power_ref = sum(abs(ho).^2);
    %epsilon = 10^(target_nmse_db / 10) * power_ref;
    

    power_ref = mean(abs(ho).^2);
    n_dirs = size(ho,1);
    epsilon = n_dirs * 10^(target_nmse_db / 10) * power_ref;



    % Solve using CVX
    cvx_begin quiet
        variable c(M) complex
        minimize( sum_square_abs(Vs' * c - hs) )
        subject to
            sum_square_abs(Vo' * c - ho) <= epsilon;

        
    cvx_end

    % if sum_square_abs(Vo' * c - ho) > epsilon
    %     warning('Constraint violated at freq bin %d: NMSE = %.2f dB', f, ...
    %         10*log10(sum_square_abs(Vo' * c - ho) / power_ref));
    %     mse = sum(abs(Vo' * c - ho).^2);
    %     power = sum(abs(ho).^2);
    %     fprintf("Actual NMSE: %.2f dB\n", 10*log10(mse / power));
    % end
    C(:, f) = c;
end
end
