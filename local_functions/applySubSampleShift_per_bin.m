function V_k_shifted = applySubSampleShift_per_bin(V_k, tau_TOA, nFFT,bin_idx)
%% applySubSampleShift
% Applies the same sub-sample time shift to all directions and channels.
%
% Inputs:
%   V_k        : [nDir x nBins x nChan] steering vectors
%   tau_TOA    : (scalar) shift in seconds (e.g., 1.5 / fs for 1.5 samples)
%   fs         : Sampling frequency in Hz
%   nFFT       : Total FFT length (BSMobj.filt_samp)
%
% Output:
%   V_shifted  : [nDir x nBins x nChan] phase-shifted vectors

    nBins = size(V_k,2);
    bin_indices = (bin_idx:nBins); % 0 for DC, 1 for fs/NFFT...
    
    % Calculate the row vector of phase shifts
    phasor = exp(-1j * 2 * pi * bin_indices * (tau_TOA / nFFT));
    
    % Apply using broadcasting
    % MATLAB will automatically apply the [1 x nBins] phasor 
    % across the [directions] and [channels] dimensions.
    V_k_shifted = V_k;
    V_k_shifted(:,bin_idx:nBins,:) = V_k(:,bin_idx:nBins,:) .* phasor;


    % [nDir, nBins, nChan] = size(V_k);
    % 
    % % 1. Create the frequency vector for the positive bins
    % % f = [0, fs/nFFT, 2*fs/nFFT, ...]
    % f = (0:nBins-1) * (fs / nFFT); % [1 x nBins]
    % 
    % % 2. Calculate the phase ramp for this specific tau
    % % Phi = -2 * pi * f * tau
    % phi = -2 * pi * f * tau_TOA; % [1 x nBins]
    % 
    % % 3. Generate the complex phase shift
    % phase_shift = exp(1j * phi); % [1 x nBins]
    % 
    % % 4. Apply to V_k
    % % V_k is [nDir x nBins x nChan]
    % % phase_shift is [1 x nBins]
    % % MATLAB handles the broadcasting across nDir and nChan automatically
    % V_shifted = V_k .* phase_shift;

end