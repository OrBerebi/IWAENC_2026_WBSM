function [c_BSM_l, c_BSM_r] = GenerateBSMfilters_v2(BSMobj, V_k, hobj_freq_grid, W)
% GenerateBSMfilters_v2.m
% Generates Binaural Synthesis Model (BSM) filters with improved structure.
%
% This version includes:
% 1. Cleaner, refactored code with helper functions for clarity.
% 2. Conditional processing for Magnitude Least Squares (MagLS).
% 3. A smooth, linear transition (cross-fade) over half an octave between
%    the Complex LS and Magnitude LS filter solutions.
%
% Inputs:
%   BSMobj          : (object) MATLAB object containing all BSM parameters.
%   V_k             : (n_mic x Q x n_freqs) Steering vectors (ATF).
%   hobj_freq_grid  : (object) HRTF data interpolated to the desired frequencies.
%   W               : (matrix) Weighting matrix for the LS solution.
%
% Outputs:
%   c_BSM_l         : (n_mic x n_freqs) BSM filters for the left ear (freq. domain).
%   c_BSM_r         : (n_mic x n_freqs) BSM filters for the right ear (freq. domain).
%
% Author: Or Berebi
% Date: 2025-10-18

%% 1. Extract Parameters from BSMobj
% This makes the main function body less cluttered.
freqs_sig = BSMobj.freqs_sig;
n_mic = BSMobj.n_mic;
normSV = BSMobj.normSV;
magLS = BSMobj.magLS;
f_cut_magLS = BSMobj.f_cut_magLS;

%% 2. Initialize Output Arrays
n_freqs = length(freqs_sig);
c_BSM_l = zeros(n_mic, n_freqs);
c_BSM_r = zeros(n_mic, n_freqs);

%% 3. Normalize Steering Vectors (if required)
% This operation is done once before the main loop for efficiency.
if normSV
    V_k = V_k ./ vecnorm(V_k, 2, 1);
end

%% 4. Main Processing Loop
% The logic is split based on whether MagLS processing is enabled.

if ~magLS
    % --- BRANCH 1: Complex LS Processing Only ---
    % If magLS is disabled, apply complex LS across all frequencies.
    fprintf('Running Complex LS processing only.\n');
    
    for f = 1:n_freqs
        V_k_curr = V_k(:, :, f);
        h_l = hobj_freq_grid.data(:, f, 1);
        h_r = hobj_freq_grid.data(:, f, 2);
        
        if any(isnan(V_k_curr), 'all')
            c_BSM_l(:, f) = zeros(n_mic, 1);
            c_BSM_r(:, f) = zeros(n_mic, 1);
        else
            c_BSM_l(:, f) = calculate_ls_solution(BSMobj, V_k_curr, h_l, W);
            c_BSM_r(:, f) = calculate_ls_solution(BSMobj, V_k_curr, h_r, W);
        end
    end
    
else
    % --- BRANCH 2: LS / MagLS with Smooth Transition ---
    % Define the half-octave transition band centered at f_cut_magLS.
    % A quarter octave corresponds to a frequency ratio of 2^0.25.
    f_ratio = 2^0.25;
    f_low = f_cut_magLS / f_ratio;
    f_high = f_cut_magLS * f_ratio;
    
    fprintf('Running LS/MagLS with smooth transition between %.1f Hz and %.1f Hz.\n', f_low, f_high);
    
    for f = 1:n_freqs
        current_freq = freqs_sig(f);
        
        % Get current frequency-dependent data
        V_k_curr = V_k(:, :, f);
        h_l = hobj_freq_grid.data(:, f, 1);
        h_r = hobj_freq_grid.data(:, f, 2);
        
        % Skip if steering vector contains NaNs
        if any(isnan(V_k_curr), 'all')
            c_BSM_l(:, f) = zeros(n_mic, 1);
            c_BSM_r(:, f) = zeros(n_mic, 1);
            continue;
        end
        
        % Get previous filter solution for MagLS initialization.
        if f == 1
            c_prev_l = zeros(n_mic, 1);
            c_prev_r = zeros(n_mic, 1);
        else
            c_prev_l = c_BSM_l(:, f-1);
            c_prev_r = c_BSM_r(:, f-1);
        end
        
        % Determine processing based on frequency band
        if current_freq < f_low
            % --- Below transition band: Pure Complex LS ---
            c_l = calculate_ls_solution(BSMobj, V_k_curr, h_l, W);
            c_r = calculate_ls_solution(BSMobj, V_k_curr, h_r, W);
            
        elseif current_freq > f_high
            % --- Above transition band: Pure Magnitude LS ---
            c_l = calculate_magls_solution(BSMobj, V_k_curr, h_l, W, c_prev_l);
            c_r = calculate_magls_solution(BSMobj, V_k_curr, h_r, W, c_prev_r);
            
        else
            % --- Inside transition band: Linear Interpolation ---
            
            % Calculate both LS and MagLS solutions
            c_ls_l = calculate_ls_solution(BSMobj, V_k_curr, h_l, W);
            c_ls_r = calculate_ls_solution(BSMobj, V_k_curr, h_r, W);
            
            c_magls_l = calculate_magls_solution(BSMobj, V_k_curr, h_l, W, c_prev_l);
            c_magls_r = calculate_magls_solution(BSMobj, V_k_curr, h_r, W, c_prev_r);
            
            % Calculate weighting factor (alpha) for cross-fading
            alpha = (current_freq - f_low) / (f_high - f_low);
            
            % Linearly interpolate between the two solutions
            c_l = (1 - alpha) * c_ls_l + alpha * c_magls_l;
            c_r = (1 - alpha) * c_ls_r + alpha * c_magls_r;
        end
        
        % Store results
        c_BSM_l(:, f) = c_l;
        c_BSM_r(:, f) = c_r;
    end
end

% c_BSM_l = conj(c_BSM_l);
% c_BSM_r = conj(c_BSM_r);
end

%% --- Helper Function for Complex Least Squares ---
function c = calculate_ls_solution(BSMobj, V_k, h, W)
% Encapsulates all logic for the Complex LS solution.

% Extract necessary parameters
Tik = BSMobj.Tik;
Jan = BSMobj.Jan;
SNR_lin = BSMobj.SNR_lin;
lambda = BSMobj.beta.^2;
gamma = 1e-3;
omega = [BSMobj.th_BSMgrid_vec, BSMobj.ph_BSMgrid_vec, BSMobj.r_BSMgrid_vec];
sector_ang = BSMobj.sector_ang;

if Tik
    c = BSM_toolbox.LeastSqueres_Tik_BSM_solution(V_k, h, lambda, gamma, omega, sector_ang);
elseif Jan
    c = BSM_toolbox.LeastSqueres_JBSM_solution(V_k, h, W);
else
    % This handles both the weighted and unweighted (Tikhonov) cases
    c = BSM_toolbox.LeastSqueres_BSM_solution(V_k, h, W, (1 / SNR_lin));
    %c = BSM_toolbox.TikhonovReg(V_k, h, (1 / SNR_lin),2);
end
end

%% --- Helper Function for Magnitude Least Squares ---
function c = calculate_magls_solution(BSMobj, V_k, h, W, c_prev)
% Encapsulates all logic for the Magnitude LS solution.

% Extract necessary parameters
Tik = BSMobj.Tik;
Jan = BSMobj.Jan;
magLS_cvx = BSMobj.magLS_cvx;
SNR_lin = BSMobj.SNR_lin;
inv_opt = BSMobj.inv_opt;
lambda = BSMobj.beta.^2;
gamma = 1e-3;
omega = [BSMobj.th_BSMgrid_vec, BSMobj.ph_BSMgrid_vec, BSMobj.r_BSMgrid_vec];
sector_ang = BSMobj.sector_ang;

if magLS_cvx
    % SDP with CVX toolbox
    c = BSM_toolbox.TikhonovReg_MagLS_CVX(V_k, h, (1 / SNR_lin), inv_opt);
else
    % Variable exchange methods
    if Tik
        c = BSM_toolbox.MagLS_Tik_BSM_solution(V_k, h, lambda, gamma, omega, sector_ang, c_prev);
    elseif Jan
        c = BSM_toolbox.MagLeastSqueres_JBSM_solution(V_k, h, c_prev, W);
    else
        % This handles both the weighted and unweighted (Tikhonov) cases
        c = BSM_toolbox.MagLeastSqueres_BSM_solution(V_k, h, c_prev, W, (1 / SNR_lin));
        
    end
end

end
