function Hnm = get_MagLS_coeff(hobj,nfft, N, Y_N, Y_high, cutOffFreq,IACC_correction)
    nfft_new = 2^11;
    hobj.shutUp = true;
    fs = hobj.fs;
    [H,zero_phase_factor]  = zero_phase_it_ambi(hobj.data,nfft_new,fs);
    
    f_vec = linspace(0,fs/2,size(H,2)).';
    Hnm = [];
    [H_l_nm_MagLS, H_r_nm_MagLS] = computeMagLS_imp(H, f_vec, N, Y_N.', cutOffFreq, Y_high.',IACC_correction);
    
    H_l_nm_MagLS = sh_freq_complement(H_l_nm_MagLS, nfft_new, 1, 2);
    H_r_nm_MagLS = sh_freq_complement(H_r_nm_MagLS, nfft_new, 1, 2);

    n_shift = 100;
    [N_fft, ch] = size(H_l_nm_MagLS);
    k = (0:N_fft-1)';  % frequency bin indices
    shift_factor = exp(-1j * 2 * pi * k * n_shift / N_fft);  % column vector

    %shift_factor = exp(-1j * 2 * pi * f_vec.' * n_shift / fs);  % [freqs x 1]

    % Apply the phase shift to each channel
    H_l_nm_MagLS = H_l_nm_MagLS .* shift_factor;  % automatically broadcast across channels
    H_r_nm_MagLS = H_r_nm_MagLS .* shift_factor;  % automatically broadcast across channels

    H_l_nm_MagLS = H_l_nm_MagLS(1:N_fft/2+1,:);
    H_r_nm_MagLS = H_r_nm_MagLS(1:N_fft/2+1,:);

    Hnm = cat(3,Hnm,H_l_nm_MagLS);
    Hnm = cat(3,Hnm,H_r_nm_MagLS);
    Hnm = permute(Hnm,[2,1,3]);

    [~,H_mls_t]   = get_p_signals(Hnm,Y_N,nfft_new,false);

    H_mls_f = zeros([size(H_mls_t,1),nfft/2+1,2]);
    H_mls_t = H_mls_t(:,1:nfft,:);

    for ear_idx = 1:2
        tmp = H_mls_t(:,:,ear_idx);
        tmp = fft(tmp,nfft,2);
        tmp = tmp(:,1:end/2+1);
        H_mls_f(:,:,ear_idx) = tmp;
    end
    Yp = pinv(Y_N);
    foo = zeros([size(Hnm,1),nfft/2+1,2]);
    for ear_idx = 1:2
        tmp = H_mls_f(:,:,ear_idx);
        foo(:,:,ear_idx) = Yp*double(tmp);
    end
    Hnm = foo;
    
   
end


