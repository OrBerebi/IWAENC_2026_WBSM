function [binaural_sig, sig_direct, sig_ambient,c_COM] = compass_bsm_static(mic_signals, Rx_global, ATF_direct, HRTF_direct, c_BSM)
%COMPASS_BSM_STATIC_DEBUG 

    %% 1. Setup
    c_BSM = c_BSM(:, [2, 1], :); % Swap Column 1 (Left) and Column 2 (Right)
    [nSamples, Q] = size(mic_signals);
    [Q_v, nBins, K] = size(ATF_direct);
    nfft = (nBins - 1) * 2;
    
    %% 2. Design Filters (Split into Direct / Ambient)
    fprintf('Designing Split Filters (N=%d)...\n', nfft);
    
    % Direct Filters (Left/Right)
    H_dir_L = zeros(nBins, Q); H_dir_R = zeros(nBins, Q);
    
    % Ambient Filters (Left/Right)
    H_amb_L = zeros(nBins, Q); H_amb_R = zeros(nBins, Q);
    
    for k = 1:nBins
        % Inputs
        V_d = squeeze(ATF_direct(:, k, :));
        h_d_L = squeeze(HRTF_direct(1, k, :)); 
        h_d_R = squeeze(HRTF_direct(2, k, :)); 
        c_BSM_L = squeeze(c_BSM(:, 1, k));      
        c_BSM_R = squeeze(c_BSM(:, 2, k));      
        Rx = Rx_global(:, :, k);
        
        % Beamformer (Wd)
        beta = 0.05 * trace(Rx) + 1e-8;
        Rx_inv = inv(Rx + beta * eye(Q));
        denom = V_d' * Rx_inv * V_d;
        
        if rcond(denom) < 1e-6
             W_d = (V_d') ./ sum(abs(V_d).^2, 1)';
        else
             W_d = denom \ (V_d' * Rx_inv);
        end
        if norm(W_d) > 5, W_d = W_d / norm(W_d) * 5; end
        
        % Residual (Wr)
        W_r = eye(Q) - V_d * W_d;
        
        % --- SEPARATE THE TERMS ---
        
        % 1. Direct Part: c_dir = Wd' * conj(h)
        % Note: h is [K x 1]. Wd is [K x Q]. Wd' is [Q x K].
        % Result is [Q x 1]
        c_d_L_k = W_d' * conj(h_d_L); 
        c_d_R_k = W_d' * conj(h_d_R);
        
        % 2. Ambient Part: c_amb = Wr' * c_BSM
        c_a_L_k = W_r' * c_BSM_L;
        c_a_R_k = W_r' * c_BSM_R;
        
        % Store
        H_dir_L(k, :) = c_d_L_k.'; H_dir_R(k, :) = c_d_R_k.';
        H_amb_L(k, :) = c_a_L_k.'; H_amb_R(k, :) = c_a_R_k.';
    end
    
    %% 3. Convert to Time Domain
    
    % Helper to convert [nBins x Q] -> [nfft x Q] impulse response
    %to_time = @(H) fftshift(real(ifft([H; conj(H(end-1:-1:2, :))], nfft, 1)), 1);
    to_time = @(H) real(ifft([H; conj(H(end-1:-1:2, :))], nfft, 1));
    
    h_dir_L_t = to_time(H_dir_L); h_dir_R_t = to_time(H_dir_R);
    h_amb_L_t = to_time(H_amb_L); h_amb_R_t = to_time(H_amb_R);
    

    tmp_l   = H_dir_L + H_amb_L; tmp_r = H_dir_R + H_amb_R;
    c_COM.f = cat(3,tmp_l,tmp_r);

    tmp_l   = h_dir_L_t + h_amb_L_t; tmp_r = h_dir_R_t + h_amb_R_t;
    c_COM.t = cat(3,tmp_l,tmp_r);
    %% 4. Apply Filters
    fprintf('Applying Filters...\n');
    
    % Apply Direct Filters
    out_d_L = zeros(nSamples, 1); out_d_R = zeros(nSamples, 1);
    for q = 1:Q
        out_d_L = out_d_L + fftfilt(h_dir_L_t(:, q), mic_signals(:, q));
        out_d_R = out_d_R + fftfilt(h_dir_R_t(:, q), mic_signals(:, q));
    end
    sig_direct = [out_d_L, out_d_R];
    
    % Apply Ambient Filters
    out_a_L = zeros(nSamples, 1); out_a_R = zeros(nSamples, 1);
    for q = 1:Q
        out_a_L = out_a_L + fftfilt(h_amb_L_t(:, q), mic_signals(:, q));
        out_a_R = out_a_R + fftfilt(h_amb_R_t(:, q), mic_signals(:, q));
    end
    sig_ambient = [out_a_L, out_a_R];
    
    % Total
    binaural_sig = sig_direct + sig_ambient;
    %binaural_sig = circshift(binaural_sig,49,1);
    
    % Normalize
    max_val = max(abs(binaural_sig(:)));
    if max_val > 0, binaural_sig = binaural_sig / max_val; end
end