function measure_atf_physical_size(ATF_single_dir, fs)
    % ATF_single_dir: [Q_mics x nBins] (The processed ATF for ONE specific direction)
    % fs: Sampling rate (e.g., 48000 or 16000)
    
    [Q, nBins] = size(ATF_single_dir);
    nfft = (nBins - 1) * 2;
    
    % 1. Reconstruct full spectrum and go to Time Domain
    ATF_full = [ATF_single_dir, conj(ATF_single_dir(:, end-1:-1:2))];
    h_time = real(ifft(ATF_full, nfft, 2));
    
    % 2. Find the arrival time (peak) for each microphone
    arrival_samples = zeros(1, Q);
    for q = 1:Q
        [~, peak_idx] = max(abs(h_time(q, :)));
        arrival_samples(q) = peak_idx;
    end
    
    % 3. Calculate Physical Metrics
    c = 343; % Speed of sound (m/s)
    
    % Max delay between any two mics in the array
    max_delay_samples = max(arrival_samples) - min(arrival_samples);
    max_delay_seconds = max_delay_samples / fs;
    max_distance_meters = max_delay_seconds * c;
    
    fprintf('\n=== ATF PHYSICAL DIMENSION TEST ===\n');
    fprintf('Max Time Delay: %.2f ms\n', max_delay_seconds * 1000);
    fprintf('Max Array Dimension: %.2f cm\n', max_distance_meters * 100);
    fprintf('Effective Array Radius: %.2f cm\n', (max_distance_meters / 2) * 100);
    fprintf('===================================\n');
    
    % Plot the impulses to visually verify
    figure('Name', 'ATF Impulse Responses', 'Color', 'w');
    plot(h_time.', 'LineWidth', 1.5);
    xlim([min(arrival_samples)-10, max(arrival_samples)+10]);
    title('Time-of-Arrival for Processed ATF');
    xlabel('Samples'); ylabel('Amplitude');
    grid on;
end