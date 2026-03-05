function p_binaural_t = bsm_binaural_renderer(p_array,c_bsm)

nfft          = max(size(c_bsm,2),size(p_array,1));
nfft          = 2^nextpow2(nfft);
c_bsm_f       = fft(c_bsm,nfft,2);
p_array_f     = fft(p_array,nfft,1);
p_binaural_l  = sum(conj(c_bsm_f(:,:,1).') .* p_array_f, 2);
p_binaural_r  = sum(conj(c_bsm_f(:,:,2).') .* p_array_f, 2);
p_binaural    = [p_binaural_l,p_binaural_r];
p_binaural_t  = ifft(p_binaural, nfft, 1, 'symmetric');

end