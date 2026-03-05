function [c_BSM_l, c_BSM_r] = GenerateBSMfilters_v3(BSMobj, V_k, hobj_freq_grid, W, toa_samples)
% GenerateBSMfilters_v2.m
% Generates Binaural Synthesis Model (BSM) filters with MagLS Phase Correction.
%
% Updates:
% - Added 'toa_samples' input to re-inject Time-of-Arrival delay into MagLS bins.
% - Prevents phase drift by separating the "Solver State" from "Output".

    %% 1. Extract Parameters
    freqs_sig   = BSMobj.freqs_sig;
    n_mic       = BSMobj.n_mic;
    normSV      = BSMobj.normSV;
    magLS       = BSMobj.magLS;
    f_cut_magLS = BSMobj.f_cut_magLS;
    
    % Default ToA if not provided (User mentioned ~26 bins)
    if nargin < 5
        if isfield(BSMobj, 'toa_samples')
            toa_samples = BSMobj.toa_samples;
        else
            toa_samples = 26; 
        end
    end

    %% 2. Initialize
    n_freqs = length(freqs_sig);
    c_BSM_l = zeros(n_mic, n_freqs);
    c_BSM_r = zeros(n_mic, n_freqs);
    
    % Determine FFT size for correct phase calculation
    % (Assuming n_freqs = nFFT/2 + 1)
    nFFT = (n_freqs - 1) * 2;

    %% 3. Normalize SV
    if normSV
        V_k = V_k ./ vecnorm(V_k, 2, 1);
    end

    %% 4. Main Processing Loop
    if ~magLS
        % --- Pure Complex LS ---
        fprintf('Running Complex LS processing only.\n');
        for f = 1:n_freqs
            V_k_curr = V_k(:, :, f);
            % Handle NaNs
            if any(isnan(V_k_curr), 'all')
                c_BSM_l(:, f) = zeros(n_mic, 1);
                c_BSM_r(:, f) = zeros(n_mic, 1);
                continue;
            end
            
            h_l = hobj_freq_grid.data(:, f, 1);
            h_r = hobj_freq_grid.data(:, f, 2);
            c_BSM_l(:, f) = calculate_ls_solution(BSMobj, V_k_curr, h_l, W);
            c_BSM_r(:, f) = calculate_ls_solution(BSMobj, V_k_curr, h_r, W);
        end
        
    else
        % --- LS / MagLS with Phase Injection ---
        f_ratio = 2^0.25;
        f_low = f_cut_magLS / f_ratio;
        f_high = f_cut_magLS * f_ratio;
        
        fprintf('Running MagLS with ToA Injection (%d samples) > %.1f Hz.\n', toa_samples, f_high);

        % Initialize "Raw" state tracking (Unshifted)
        c_prev_l_raw = zeros(n_mic, 1);
        c_prev_r_raw = zeros(n_mic, 1);

        for f = 1:n_freqs
            current_freq = freqs_sig(f);
            V_k_curr = V_k(:, :, f);
            
            if any(isnan(V_k_curr), 'all')
                c_BSM_l(:, f) = zeros(n_mic, 1);
                c_BSM_r(:, f) = zeros(n_mic, 1);
                c_prev_l_raw = zeros(n_mic, 1); % Reset state
                c_prev_r_raw = zeros(n_mic, 1);
                continue;
            end

            h_l = hobj_freq_grid.data(:, f, 1);
            h_r = hobj_freq_grid.data(:, f, 2);

            % --- Calculate Phase Shift Term for MagLS ---
            % exp(-j * 2pi * f * delay)
            % f is bin index (0 to N/2)
            bin_idx = f - 1; 
            phase_val = 2 * pi * bin_idx * (toa_samples / nFFT);
            phasor = exp(1j * phase_val);

            if current_freq < f_low
                % --- Low Freq: Complex LS (Has Natural Delay) ---
                c_l = calculate_ls_solution(BSMobj, V_k_curr, h_l, W);
                c_r = calculate_ls_solution(BSMobj, V_k_curr, h_r, W);
                
                % Update RAW state for next MagLS bin
                c_prev_l_raw = c_l;
                c_prev_r_raw = c_r;
                
            elseif current_freq > f_high
                % --- High Freq: MagLS (Needs Delay Injection) ---
                
                % 1. Calculate Raw MagLS (targets phase of c_prev_raw)
                c_l_raw = calculate_magls_solution(BSMobj, V_k_curr, h_l, W, c_prev_l_raw);
                c_r_raw = calculate_magls_solution(BSMobj, V_k_curr, h_r, W, c_prev_r_raw);
                
                % 2. Update State (Keep it unshifted to prevent accumulation!)
                c_prev_l_raw = c_l_raw;
                c_prev_r_raw = c_r_raw;
                
                % 3. Apply ToA Shift to OUTPUT only
                c_l = c_l_raw * phasor;
                c_r = c_r_raw * phasor;

                % c_l = c_l_raw;
                % c_r = c_r_raw;
                
            else
                % --- Transition: Crossfade ---
                
                % 1. LS Part (Has Natural Delay)
                c_ls_l = calculate_ls_solution(BSMobj, V_k_curr, h_l, W);
                c_ls_r = calculate_ls_solution(BSMobj, V_k_curr, h_r, W);
                
                % 2. MagLS Part (Raw / Unshifted)
                c_mag_l_raw = calculate_magls_solution(BSMobj, V_k_curr, h_l, W, c_prev_l_raw);
                c_mag_r_raw = calculate_magls_solution(BSMobj, V_k_curr, h_r, W, c_prev_r_raw);
                
                % 3. Calculate Crossfade Alpha
                alpha = (current_freq - f_low) / (f_high - f_low);
                
                % --- THE FIX: State Tracking ---
                % Strip the natural delay from LS to match the "Raw" domain
                c_ls_l_raw = c_ls_l * conj(phasor);
                c_ls_r_raw = c_ls_r * conj(phasor);
                
                % Update State using ONLY unshifted (Raw) signals
                c_prev_l_raw = (1 - alpha) * c_ls_l_raw + alpha * c_mag_l_raw;
                c_prev_r_raw = (1 - alpha) * c_ls_r_raw + alpha * c_mag_r_raw;
                
                % --- THE OUTPUT MIX ---
                % Apply ToA Shift to MagLS Part BEFORE Mixing for output
                c_mag_l_shifted = c_mag_l_raw * phasor;
                c_mag_r_shifted = c_mag_r_raw * phasor;
                
                % Mix (Both signals now safely have the delay)
                c_l = (1 - alpha) * c_ls_l + alpha * c_mag_l_shifted;
                c_r = (1 - alpha) * c_ls_r + alpha * c_mag_r_shifted;
            end
            
            % Store Final Output
            c_BSM_l(:, f) = c_l;
            c_BSM_r(:, f) = c_r;
        end
    end

    %% 5. Global Gain Matching (Tikhonov De-biasing)
    % fprintf('Applying global gain matching across all frequencies and ears...\n');
    % [c_BSM_l, c_BSM_r, alpha_gain] = match_bsm_magls_gain(c_BSM_l, c_BSM_r, V_k,  hobj_freq_grid.data, freqs_sig, f_cut_magLS);

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
