function [p_f,p_t] = get_p_signals(hnm,Y,nfft,flag_circ)
    if nargin<4
        flag_circ = false;
    end
    p_f = [];
    p_hat_omega_l = Y*hnm(:,:,1);  
    p_hat_omega_r = Y*hnm(:,:,2);
    p_f = cat(3,p_f,p_hat_omega_l);
    p_f = cat(3,p_f,p_hat_omega_r);
    p_f_tmp = p_f;
    p_f_tmp(:, end+1 : nfft, :) = 0;
    p_t = ifft(p_f_tmp, [], 2, 'symmetric');
    
    if flag_circ
        p_t  = circshift(p_t, floor(nfft / 2), 2);
        threshold_dB = 45;
        [p_t, min_delay, start_samples] = alignHRIRToMinimumDelay(p_t, threshold_dB);
    end

end