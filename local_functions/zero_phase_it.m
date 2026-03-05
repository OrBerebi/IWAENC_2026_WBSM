function h_t = zero_phase_it(h_ref,nfft,fs)


    
    h_ref_tmp       = permute(h_ref,[2,1,3]);
    group_delay_ref = group_delay_multi_channel(h_ref_tmp,nfft);
    % === 1. FFT and truncate to positive frequencies ===
    H_ref = fft(h_ref_tmp, nfft, 1);
    H_ref = H_ref( 1:nfft/2+1, :,:);
 

    % === 2. Frequency axis ===
    freqs = (0:nfft/2)' * (fs / nfft);  % in Hz

    % === 3. Frequency mask to limit to stable region ===
    %max_fit_freq = 1.3e3;  % Hz
    max_fit_freq = 8e3;  % Hz
    %max_fit_freq = 0.25e3;  % Hz
    fit_mask1 = freqs < max_fit_freq;
    %min_fit_freq = 0.4e3;  % Hz
    min_fit_freq = 20;  % Hz
    fit_mask2 = freqs > min_fit_freq;
    fit_mask = logical((fit_mask1) .* (fit_mask2));

    % === 4. Compute group delay difference ===
    delta_group_delay = group_delay_ref;  % [dirs x freqs x ears]

    % === 6. Estimate global delay ===
    estimated_delay = mean(mean(delta_group_delay(fit_mask),2));  % scalar in samples
    %estimated_delay = 65;
    fprintf('Estimated relative group delay (Left): %.2f samples\n', estimated_delay);

    % === 7. Apply linear phase correction ===
    phase_correction = exp(1j * 2 * pi * freqs * estimated_delay / fs);  % [freqs x 1]

    % Apply correction to H_ref
    H_ref = H_ref .* phase_correction;
   

    % Preallocate full Hermitian spectrum
    H_full = zeros(nfft, size(H_ref,2),size(H_ref,3));

    % Copy positive frequencies
    H_full(1:nfft/2+1, :,:) = H_ref;

    % Mirror to negative frequencies: exclude DC and Nyquist
    H_full(nfft/2+2:end, :,:) = conj(flip(H_ref(2:nfft/2, :,:), 1));

    % Perform real-valued IFFT using Hermitian symmetry
    h_t = ifft(H_full, nfft, 1, 'symmetric');
    h_t = permute(h_t,[2,1,3]);



    

    % norm to max
    %h_t = h_t./max(max(abs(h_t)));

   



end