function Rx = estimate_Rx(mic_sig,stft_params,smoothing)
% This function calculate the Signal Covariance Matrix SCM using STFT on
% the time domain signal. 

%   Initialize parameters
nfft       = stft_params.nfft;
hop        = stft_params.hop;
%window     = bartlett(nfft);
window     = stft_params.win;
[~, M]     = size(mic_sig);
desired_fs = stft_params.fs;
nFreq      = nfft/2 + 1;
pt         = mic_sig;

% Perform STFT
% Explicitly label the sampling frequency
%[S, F, T] = stft(pt, 'Window', window, 'OverlapLength', (window - hop), 'FFTLength', nfft, 'Fs', desired_fs);
[S, ~, ~] = stft(pt, desired_fs, 'Window', window, 'OverlapLength', (length(window) - hop), 'FFTLength', nfft);
%[S, ~, ~] = stft(pt, window, hop, nfft, desired_fs);
p_STFT    = S(1:size(S,1) / 2 + 1,:,:); % discard negative frequencies, +1 is the DC component
[~,T,~]   = size(p_STFT);
p_STFT    = permute(p_STFT,[3 2 1]);

% initialize arrays
Rx  = zeros([M M nFreq]);
K   = 5;
win = bartlett(2 * K + 1);

% SCM Loop
if smoothing
    for f=1:nFreq % Each freq bin
        for t=1:T  % Each time bin
            temp  = zeros([M M]);
            sumup = 0;
            for k=-K:K
                if (f+k>0 && f+k<(nFreq+1)) 
                    temp  = temp + win(k + K + 1) * p_STFT(:,t,f + k) * p_STFT(:,t,f + k)';
                    sumup = sumup + win(k + K + 1);
                end 
            end
            temp      = (1/sumup)*temp;
            Rx(:,:,f) = Rx(:,:,f) + (1/T) * temp;
%             Rx(:,:,f)=V_k(:,:,f)*V_k(:,:,f)';     % Case Rs = I*ones
        end
    end

% No Smoothing (Time domain)
else
    for f=1:nFreq
        Rx(:,:,f) = (1/T) * p_STFT(:,:,f) * p_STFT(:,:,f)';
    end
end

end
