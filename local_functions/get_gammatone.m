function [BSD_f,f_c] = get_gammatone(BRIR_f_ref,BRIR_f_test,f_band,fs,ERBstruct)
    % BRIR_f_ref: [time_samples x ears (2)]
    % BRIR_f_test: [time_samples x ears (2)]
    % f_band = [20, 20e3] example for the freqncy range of the filtering
    % fs = 48e3

    BRIR_f_ref = permute(BRIR_f_ref,[2,3,1]);
    BRIR_f_ref_L = squeeze(BRIR_f_ref(:,1,:));
    BRIR_f_ref_R = squeeze(BRIR_f_ref(:,2,:));
    BRIR_f_test = permute(BRIR_f_test,[2,3,1]);
    BRIR_f_test_L = squeeze(BRIR_f_test(:,1,:));
    BRIR_f_test_R = squeeze(BRIR_f_test(:,2,:));
    
    nfft = (size(BRIR_f_ref_L,1)-1)*2;
    
    % transform to time
    BRIR_f_ref_L(end+1:nfft,:) = 0;
    BRIR_f_test_L(end+1:nfft,:) = 0;
    BRIR_t_ref_L = ifft(BRIR_f_ref_L, nfft, 1, "symmetric");
    BRIR_t_test_L = ifft(BRIR_f_test_L, nfft, 1, "symmetric");

    BSD_f_L = BSD_ERB_fast(BRIR_t_ref_L, BRIR_t_test_L, ERBstruct);


    BRIR_f_ref_R(end+1:nfft,:) = 0;
    BRIR_f_test_R(end+1:nfft,:) = 0;
    BRIR_t_ref_R = ifft(BRIR_f_ref_R, nfft, 1, "symmetric");
    BRIR_t_test_R = ifft(BRIR_f_test_R, nfft, 1, "symmetric");

    BSD_f_R = BSD_ERB_fast(BRIR_t_ref_R, BRIR_t_test_R, ERBstruct);

    f_c = ERBstruct.f_c;
    BSD_f = cat(3,BSD_f_L,BSD_f_R);

end