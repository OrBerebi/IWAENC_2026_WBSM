function [H_ref,phase_correction] = zero_phase_it_ambi(h_ref,nfft,fs)

    group_delay_ref = group_delay_multi_channel(h_ref(:,:,1).',nfft);
    % === 1. FFT and truncate to positive frequencies ===
    H_ref = fft(h_ref, nfft, 2);
    H_ref = H_ref(:, 1:nfft/2+1, :);
 

    % === 2. Frequency axis ===
    freqs = (0:nfft/2)' * (fs / nfft);  % in Hz

    % === 3. Frequency mask to limit to stable region ===
    max_fit_freq = 0.25e3;  % Hz
    fit_mask = freqs < max_fit_freq;

    % === 4. Compute group delay difference ===
    delta_group_delay = group_delay_ref;  % [dirs x freqs x ears]

    % === 5. Average over directions and ears ===
    delta_group_delay_avg = squeeze(mean(delta_group_delay, 2));  % [freqs x 1]

    % === 6. Estimate global delay ===
    estimated_delay = mean(delta_group_delay_avg(fit_mask));  % scalar in samples

    fprintf('Estimated relative group delay: %.2f samples\n', estimated_delay);

    % === 7. Apply linear phase correction ===
    phase_correction = exp(1j * 2 * pi * freqs * estimated_delay / fs);  % [freqs x 1]

    % Apply correction to H_ref
    for ear = 1:2
        for dir = 1:size(H_ref, 1)
            H_ref(dir,:,ear) = squeeze(H_ref(dir,:,ear)) .* phase_correction.';
        end
    end

end