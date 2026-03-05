function c_BSM = get_BSM_filters(BSMobj,hobj_BSM,V_k,W)



nFFT = (size(V_k,3)-1)*2;
nfft_new = nFFT;

tmp = V_k;
tmp(:,:,end+1 : nfft_new) = 0;
V_t = ifft(tmp, [], 3, 'symmetric');
V_t(:,:,nfft_new/2:end) = 0;
group_delays_L = group_delay_multi_channel(squeeze(V_t(2,:,:)).',nFFT);
[~,f_idx_min] = min(abs(BSMobj.freqs_sig - 2e3));
tau_TOA_V= -1*mean(mean(group_delays_L(1:f_idx_min,:),1));
%V_k = fft(V_t,nFFT,3);
%V_k = V_k(:,:,1:nFFT/2+1);

% for f = 1:size(V_k,3)
%     bin_idx = f - 1; 
%     phase_val = 2 * pi * bin_idx * (tau_TOA_V / nFFT);
%     phasor = exp(-1j * phase_val);
%     V_k = V_k * phasor;
% end


hobj = hobj_BSM;
hobj = hobj.toTime(nFFT,true);
[directions, time_samples, ears] = size(hobj.data);
pad_length  = nfft_new - time_samples;
hrir_padded = cat(2, hobj.data, zeros(directions, pad_length, ears));

group_delays_L = group_delay_multi_channel(hrir_padded(:,:,1).',nFFT);

[~,f_idx_min] = min(abs(BSMobj.freqs_sig - 2e3));
tau_TOA= -1*mean(mean(group_delays_L(1:f_idx_min,:),1));

if BSMobj.arrayType == 6
    tau_TOA = tau_TOA - tau_TOA_V;
else
    tau_TOA = tau_TOA;
end
%tau_TOA = 0;

% figure
% f_vec = linspace(0,hobj.fs/2,nFFT/2+1).';
% semilogx(f_vec,group_delays_L)
% xlim([100,5e3]);
% xlabel("freq")
% ylabel("samples")

H = hrir_padded;
%H = circshift(hrir_padded, -26, 2);



H = fft(H, nfft_new, 2);
H = H(:, 1:nfft_new/2+1, :);

%H = applySubSampleShift(H, tau_TOA, BSMobj.fs, nFFT);

hobj_BSM.data = H;

n_shift = nfft_new/2;
toa_samples = -1*(tau_TOA); % Average ToA delay to preserve phase structure
%toa_samples = 0; % Average ToA delay to preserve phase structure

%n_shift = 26;

% Complex version
BSMobj.magLS = false;
BSMobj.magLS_cvx = false;
tic;

[c_BSM_cmplx_l, c_BSM_cmplx_r] = BSM_toolbox.GenerateBSMfilters_v3(BSMobj, V_k, hobj_BSM,W,toa_samples);
c_BSM_cmplx_l = phase_shift_BSM(c_BSM_cmplx_l, nfft_new, n_shift);
c_BSM_cmplx_r = phase_shift_BSM(c_BSM_cmplx_r, nfft_new, n_shift);
c_BSM_f_cmplx = cat(3,c_BSM_cmplx_l,c_BSM_cmplx_r);

%%======Post-processing BSM filters (time domain)
[c_BSM_cmplx_l_time_cs, c_BSM_cmplx_r_time_cs] = ...
    BSM_toolbox.PostProcessBSMfilters(BSMobj, c_BSM_cmplx_l, c_BSM_cmplx_r);

c_BSM_t_cmplx = cat(3, c_BSM_cmplx_l_time_cs, c_BSM_cmplx_r_time_cs);

elapsed_time_LS = toc;
%fprintf('LS execution time: %.6f seconds\n', elapsed_time_LS);





% MagLS version
BSMobj.magLS = true;
BSMobj.magLS_cvx = false;
tic;
[c_BSM_mag_l, c_BSM_mag_r] = BSM_toolbox.GenerateBSMfilters_v3(BSMobj, V_k, hobj_BSM,W,toa_samples);

% for f_idx = 1:nFFT/2 +1 
%     c = c_BSM_mag_l(:,f_idx);
%     V = V_k(:,:,f_idx);
% 
%     h_hat(:,f_idx) = squeeze(c'*V);
% end
% h_hat(:,end+1 : nFFT) = 0;
% h_hat_t = ifft(h_hat, [], 2, 'symmetric');
% bsm_group_delays_L = group_delay_multi_channel(h_hat_t(:,:).',nFFT);
% figure
% f_vec = linspace(0,hobj.fs/2,nFFT/2+1).';
% semilogx(f_vec,bsm_group_delays_L)
% xlim([100,5e3]);
% xlabel("freq")
% ylabel("samples")


c_BSM_mag_l = phase_shift_BSM(c_BSM_mag_l, nfft_new, n_shift);
c_BSM_mag_r = phase_shift_BSM(c_BSM_mag_r, nfft_new, n_shift);
c_BSM_f_mls   = cat(3,c_BSM_mag_l,c_BSM_mag_r);


[c_BSM_mag_l_time_cs, c_BSM_mag_r_time_cs] = ...
    BSM_toolbox.PostProcessBSMfilters(BSMobj, c_BSM_mag_l, c_BSM_mag_r);
c_BSM_t_mag   = cat(3, c_BSM_mag_l_time_cs, c_BSM_mag_r_time_cs);

elapsed_time_MagLS = toc;
%fprintf('MagLS execution time: %.6f seconds\n', elapsed_time_MagLS);



c_BSM.t_ls    = c_BSM_t_cmplx;
c_BSM.t_mls   = c_BSM_t_mag;
c_BSM.f_ls    = c_BSM_f_cmplx;
c_BSM.f_mls   = c_BSM_f_mls;


%c_BSM.f_ls = apply_covariance_constraint(c_BSM.f_ls, V_k, hobj_BSM.data);
%c_BSM.f_mls = apply_covariance_constraint(c_BSM.f_mls, V_k, hobj_BSM.data);

%c_BSM.f_ls   = applySubSampleShift(c_BSM.f_ls , tau_TOA,  nFFT);

%c_BSM.f_mls  = applySubSampleShift(c_BSM.f_mls , 2*tau_TOA, nFFT);

f_ratio = 2^0.25;
f_high = BSMobj.f_cut_magLS * f_ratio;
[~,bin_idx] = min(abs(BSMobj.freqs_sig - f_high));

%c_BSM.f_ls   = applySubSampleShift_per_bin(c_BSM.f_ls, -1*tau_TOA, nFFT,bin_idx);

%c_BSM.f_mls  = applySubSampleShift_per_bin(c_BSM.f_mls, 2*tau_TOA, nFFT,bin_idx+1);




end






function Hnm_shifted = phase_shift_BSM(Hnm_pos, nfft_new, n_shift)
% Apply a time-domain shift (in samples) to BSM SH coefficients
% 
% Hnm_pos: [nCoeff × nFreq_pos × nEars], positive frequencies only
% nfft_new: FFT length (full spectrum)
% n_shift: shift in samples
%
% Output:
%   Hnm_shifted: [nCoeff × nFreq_pos × nEars] phase-shifted coefficients

    [nCoeff, nFreq_pos] = size(Hnm_pos);

    % --- 1) Construct full spectrum (positive + negative freqs)
    % Drop the DC (1st) and Nyquist (last) when mirroring to avoid duplication
    Hnm_full = cat(2, ...
        Hnm_pos, ...
        conj(flip(Hnm_pos(:,2:end-1,:),2)) );   % [nCoeff × N_fft × nEars]

    % --- 2) Apply phase shift factor
    k = (0:nfft_new-1)';                                    % [N_fft × 1]
    shift_factor = exp(-1j * 2 * pi * k * n_shift / nfft_new);  % col vector
    shift_factor = reshape(shift_factor, [1 nfft_new 1]);   % broadcast

    Hnm_full = Hnm_full .* shift_factor;  % [nCoeff × N_fft × nEars]

    % --- 3) Keep only positive frequencies again
    Hnm_shifted = Hnm_full(:,1:nFreq_pos);

end
