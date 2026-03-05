function y_eq = equalize_loudness_rms(y, target)
% EQUALIZE_LOUDNESS_RMS  Equalize loudness (RMS) of a stereo signal.
%
% Usage:
%   y_eq = equalize_loudness_rms(y)
%   y_eq = equalize_loudness_rms(y, target)
%
% Inputs:
%   y       : stereo input signal (N x 2)
%   target  : (optional)
%             - scalar: desired RMS value
%             - stereo signal (N x 2): target signal whose RMS to match
%
% Output:
%   y_eq    : loudness-equalized stereo signal
%
% Example:
%   y_eq = equalize_loudness_rms(y, 0.1);
%   y_eq = equalize_loudness_rms(y, ref_signal);

    % Check input
    if size(y,2) ~= 2
        error('Input signal y must be stereo (N x 2).');
    end

    % Compute RMS per channel
    rms_y = sqrt(mean(y.^2, 1));

    % Determine target RMS
    if nargin < 2
        target_rms = mean(rms_y);  % keep same average loudness
    elseif isscalar(target)
        target_rms = target * [1 1];  % same for both channels
    elseif ismatrix(target) && size(target,2) == 2
        target_rms = sqrt(mean(target.^2, 1)); % match to target signal RMS
    else
        error('Target must be a scalar or stereo signal (N x 2).');
    end

    % Compute gain per channel
    gain = target_rms ./ rms_y;

    % Apply gain
    y_eq = y .* gain;

    % Optional: clip if needed
    y_eq = max(min(y_eq, 1), -1);
end
