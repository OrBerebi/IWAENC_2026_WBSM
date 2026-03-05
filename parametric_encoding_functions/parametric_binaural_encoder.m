function [binaural_sig, binaural_sig_s, binaural_sig_d, coherence_debug] = parametric_binaural_encoder(mic_signals, ATF_sources, ATF_ambient, Y_sources, Y_ambient, D_bin, stft_params)
%PARAMETRIC_BINAURAL_ENCODER Parametric (COMPASS) Binaural Encoder
%
%   FINAL ROBUST VERSION:
%   1. Separates output into Total, Source (Direct), and Ambient streams.
%   2. Uses Robust Matched Filter (Delay-and-Sum) for stability.
%   3. Includes robust frequency resampling (Mag/Phase separate).
%   4. INCLUDES COHERENCE DEBUGGING monitor and plot.
%
%   OUTPUTS:
%     binaural_sig    : [nSamples x 2] Total Mix
%     binaural_sig_s  : [nSamples x 2] Source Stream Only
%     binaural_sig_d  : [nSamples x 2] Ambient Stream Only
%     coherence_debug : [nBins x nFrames x K] Spatial Coherence Score
%
%   Reference:
%   [1] McCormack et al., "Parametric Ambisonic Encoding of Arbitrary
%       Microphone Arrays", IEEE/ACM TASLP, 2022.

    %% 1. Unpack Parameters & Setup Targets
    [nSamples, Q] = size(mic_signals);
    
    nfft = stft_params.nfft;
    hop  = stft_params.hop;
    win  = stft_params.win;
    
    % The target number of bins for processing
    nBins_target = nfft/2 + 1;
    
    %% 2. Input Validation & Independent Resampling
    
    % Get dimensions of input data
    [Q_src, nBins_src, K] = size(ATF_sources);
    [Q_amb, nBins_amb, L] = size(ATF_ambient);
    [nEars, nBins_dec, M_sh] = size(D_bin);
    
    % Consistency Check: Channel Counts
    if Q_src ~= Q || Q_amb ~= Q
        error('Mismatch: ATF rows (%d) must match mic channels (%d).', Q_src, Q);
    end
    
    % --- Helper Function for Resampling ---
    resample_data = @(M, nIn, nOut) perform_resampling(M, nIn, nOut);
    
    % --- Resample ATF_sources ---
    if nBins_src ~= nBins_target
        fprintf('Resampling Sources: %d -> %d bins...\n', nBins_src, nBins_target);
        ATF_sources = resample_data(ATF_sources, nBins_src, nBins_target);
    end
    
    % --- Resample ATF_ambient ---
    if nBins_amb ~= nBins_target
        fprintf('Resampling Ambient: %d -> %d bins...\n', nBins_amb, nBins_target);
        ATF_ambient = resample_data(ATF_ambient, nBins_amb, nBins_target);
    end
    
    % --- Resample D_bin ---
    if nBins_dec ~= nBins_target
        fprintf('Resampling Decoder: %d -> %d bins...\n', nBins_dec, nBins_target);
        D_bin = resample_data(D_bin, nBins_dec, nBins_target);
    end
    
    %% 3. Initialization
    % Pre-compute Diffuse Coherence Matrix (D) per bin
    D_diffuse = zeros(Q, Q, nBins_target);
    for f = 1:nBins_target
        A_grid_f = squeeze(ATF_ambient(:, f, :)); % Q x L
        D_diffuse(:, :, f) = (A_grid_f * A_grid_f') / L;
    end
    
    % Covariance Averaging Memory
    alpha = 0.97; 
    Cx_mem = zeros(Q, Q, nBins_target);
    
    % Output Buffers (Overlap-Add)
    nFrames = ceil((nSamples) / hop);
    padded_len = (nFrames * hop) + nfft;
    mic_signals_pad = [mic_signals; zeros(padded_len - nSamples, Q)];
    
    binaural_out_buffer   = zeros(padded_len, 2); 
    binaural_out_buffer_s = zeros(padded_len, 2); 
    binaural_out_buffer_d = zeros(padded_len, 2); 
    
    % DEBUG: Coherence Monitor
    coherence_debug = zeros(nBins_target, nFrames, K);
    
    %% 4. Processing Loop
    
    fprintf('Parametric Encoding: Processing %d frames...\n', nFrames);
    
    for frame_idx = 0:nFrames-1
        
        % --- STFT Analysis ---
        idx_start = frame_idx * hop + 1;
        idx_end = idx_start + nfft - 1;
        
        if idx_end > size(mic_signals_pad, 1), break; end
        
        x_frame = mic_signals_pad(idx_start:idx_end, :); 
        x_win = x_frame .* win; 
        X_f = fft(x_win, nfft, 1); 
        X_half = X_f(1:nBins_target, :).'; % (Q x nBins_target)
        
        Bin_frame   = zeros(2, nBins_target);
        Bin_frame_s = zeros(2, nBins_target);
        Bin_frame_d = zeros(2, nBins_target);
        
        % --- Frequency Bin Processing ---
        for k = 1:nBins_target
            
            % Load Bin Data
            x_k = X_half(:, k);                     % (Q x 1)
            As = squeeze(ATF_sources(:, k, :));     % (Q x K)
            Ad = squeeze(ATF_ambient(:, k, :));     % (Q x L)
            D_diff = squeeze(D_diffuse(:, :, k));   % (Q x Q)
            H_bin = squeeze(D_bin(:, k, :));        % (2 x M)
            
            % --- DEBUG: COHERENCE CHECK ---
            norm_x = norm(x_k);
            if norm_x > 1e-12
                % Compute spatial correlation between Model (As) and Data (x_k)
                dot_prod = abs(As' * x_k);
                norm_As = sqrt(sum(abs(As).^2, 1)).'; % (K x 1)
                coh_val = dot_prod ./ (norm_As * norm_x);
                coherence_debug(k, frame_idx+1, :) = coh_val;
            end
            
            % Update Covariance (Cx) - Kept for optional LCMP use
            cx_inst = x_k * x_k';
            Cx_mem(:, :, k) = alpha * Cx_mem(:, :, k) + (1-alpha) * cx_inst;
            Cx = Cx_mem(:, :, k); 
            
            % -------------------------------------------------------------
            % 1. Source Stream Rendering
            % -------------------------------------------------------------
            
            % OPTION A: Matched Filter (Delay-and-Sum) - ACTIVE
            % vec_energies = sum(abs(As).^2, 1).'; % (K x 1)
            % vec_energies(vec_energies < 1e-12) = 1e-6;
            % Ws = (1 ./ vec_energies) .* As'; % (K x Q)
            
            % OPTION B: LCMP (Commented Out)
            tr_Cx = trace(Cx);
            beta = 0.05 * tr_Cx + 1e-7; 
            R = Cx + beta * eye(Q);
            Ris_As = R \ As; 
            denom = As' * Ris_As;
            Ws = pinv(denom) * Ris_As'; 
            
            % Extract Sources
            s_est = Ws * x_k;
            a_s = Y_sources * s_est;
            
            % -------------------------------------------------------------
            % 2. Ambient Stream Rendering (Residual + PWD)
            % -------------------------------------------------------------
            
            Wd = eye(Q) - As * Ws;
            d_res = Wd * x_k;
         
            % Energy-Preserving SVD
            [U_d, ~, V_d] = svd(Ad', 'econ');
            if size(U_d, 2) > Q, U_d = U_d(:, 1:Q); end
            
            A_hat_d = (1/sqrt(L)) * U_d * V_d';
            z_d = A_hat_d * d_res;
            
            % Decorrelation
            rand_phase = exp(1i * (2*pi*rand(L,1) - pi));
            z_d = z_d .* rand_phase;
            
            % Equalization & Encoding
            tr_D = trace(D_diff);
            if tr_D < 1e-12, tr_D = 1e-6; end
            Ed = 1 / sqrt(tr_D);
            
            a_d = Ed * Y_ambient * z_d;
            
            % -------------------------------------------------------------
            % 3. Combination & Output
            % -------------------------------------------------------------
            
            a_par = a_s + a_d;
            
            Bin_frame(:, k)   = H_bin * a_par;
            Bin_frame_s(:, k) = H_bin * a_s;
            Bin_frame_d(:, k) = H_bin * a_d;
        end
        
        % --- ISTFT Synthesis ---
        
        % Total
        Bin_frame_full = [Bin_frame, conj(Bin_frame(:, end-1:-1:2))];
        b_time = real(ifft(Bin_frame_full, nfft, 2));
        binaural_out_buffer(idx_start:idx_end, :) = ...
            binaural_out_buffer(idx_start:idx_end, :) + b_time.';
            
        % Source
        Bin_frame_full_s = [Bin_frame_s, conj(Bin_frame_s(:, end-1:-1:2))];
        b_time_s = real(ifft(Bin_frame_full_s, nfft, 2));
        binaural_out_buffer_s(idx_start:idx_end, :) = ...
            binaural_out_buffer_s(idx_start:idx_end, :) + b_time_s.';
            
        % Ambient
        Bin_frame_full_d = [Bin_frame_d, conj(Bin_frame_d(:, end-1:-1:2))];
        b_time_d = real(ifft(Bin_frame_full_d, nfft, 2));
        binaural_out_buffer_d(idx_start:idx_end, :) = ...
            binaural_out_buffer_d(idx_start:idx_end, :) + b_time_d.';
    end
    
    % Trim Buffers
    binaural_sig   = binaural_out_buffer(1:nSamples, :);
    binaural_sig_s = binaural_out_buffer_s(1:nSamples, :);
    binaural_sig_d = binaural_out_buffer_d(1:nSamples, :);
    
    %% 5. Debug Visualization
    figure('Name', 'Source Model Coherence Check', 'Color', 'w');
    
    % Plot coherence for Source 1
    src_idx = 1;
    coh_map = squeeze(coherence_debug(:, :, src_idx));
    
    % Time and Freq Axes
    t_axis = (0:nFrames-1) * hop / 48000; % Approx seconds
    f_axis = linspace(0, 24000, nBins_target);
    
    imagesc(t_axis, f_axis, coh_map);
    axis xy; colormap jet; colorbar;
    clim([0 1]); 
    
    title(sprintf('Spatial Coherence: Source %d vs Mic Signal', src_idx));
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    
    avg_coh = mean(coh_map(:));
    fprintf('\n--- COHERENCE REPORT ---\n');
    fprintf('Global Average Coherence: %.4f\n', avg_coh);
    if avg_coh < 0.2
        fprintf('CRITICAL: Very low coherence (<0.2). Indicates MODEL MISMATCH.\n');
    else
        fprintf('STATUS: Coherence acceptable.\n');
    end
    fprintf('------------------------\n');
end

function M_out = perform_resampling(M_in, nIn, nOut)
    % ROBUST RESAMPLING: Interpolates Magnitude and Unwrapped Phase separately.
    freq_in  = linspace(0, 1, nIn);
    freq_out = linspace(0, 1, nOut);
    
    M_perm = permute(M_in, [2, 1, 3]);
    M_mag = abs(M_perm);
    M_phi = unwrap(angle(M_perm)); 
    
    M_mag_new = interp1(freq_in, M_mag, freq_out, 'linear', 'extrap');
    M_phi_new = interp1(freq_in, M_phi, freq_out, 'linear', 'extrap');
    
    M_interp = M_mag_new .* exp(1i * M_phi_new);
    M_out = permute(M_interp, [2, 1, 3]);
end