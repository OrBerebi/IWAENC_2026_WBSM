function group_delays = group_delay_multi_channel(impulse_response, N_fft)
    % group_delay_multi_channel Compute group delay for multi-channel impulse responses
    %
    % Inputs:
    %   impulse_response : [N x ch] time-domain impulse responses
    %   N_fft (optional) : number of FFT points (zero-padding if N_fft > N)
    %
    % Output:
    %   group_delays     : [N_fft/2+1 x ch] group delay per frequency bin and channel

    if nargin < 2
        N_fft = size(impulse_response, 1);  % No padding by default
    end

    [N, ch] = size(impulse_response);

    % Time index vector
    n = (0:N-1)';
    nh = impulse_response .* n;

    % Compute FFT with zero-padding to N_fft
    dft_h = fft(impulse_response, N_fft, 1) / N;
    dft_nh = fft(nh, N_fft, 1) / N;

    % Compute group delay
    eps_val = 1e-8;
    group_delays = real(dft_nh ./ (dft_h + eps_val));

    % Keep only the positive frequency bins
    group_delays = group_delays(1:floor(N_fft/2)+1, :);

    % Check for NaNs
    if any(isnan(group_delays), 'all')
        warning('NaNs detected in group delay!');
        min_val = min(abs(dft_h), [], 'all');
        fprintf('Min abs(dft_h): %g\n', min_val);
        error('NaNs in group delay computation');
    end
end