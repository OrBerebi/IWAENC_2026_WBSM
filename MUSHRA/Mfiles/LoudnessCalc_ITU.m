function loudness = LoudnessCalc_ITU(audioData)
%
%--- Loudness calculation according to ITU-R BS. 1770 ---
%
% This script calculates the loudness value for a mono or stereo .wav-file
% according to ITU-R BS. 1770
%
% Input: RIRwavFile     - name of a BRIR .wav-file on the form: 'name.wav'.
%        StimuliwavFile - name of a stimuli .wav-file on the form: 'name.wav'.
%
% Output: The calculated loudness value.
%
% WARNING: The filters coefficients in this code are for a sampling rate of 48 kHz!
%
% Zamir Ben-Hur
% 29.7.2015
%

% Filter coefficients for the modeling of the acoustic effects of the head (48kHz)
Bhead = [1.53512485958697 -2.69169618940638 1.19839281085285];
Ahead = [1 -1.69065929318241 0.73248077421585];

% RLB filter coefficients (48kHz)
BRLB = [1 -2 1];
ARLB = [1 -1.99004745483398 0.99007225036621];

% Filtering and calculation according to ITU-R BS. 1770-1 (stereo)
af1ch1 = filter(Bhead, Ahead, audioData(:,1));
af1ch2 = filter(Bhead, Ahead, audioData(:,2));

af1ch1filt = filter(BRLB, ARLB, af1ch1(:));
af1ch2filt = filter(BRLB, ARLB, af1ch2(:));

af1ch1filtsq = af1ch1filt.^2;
af1ch2filtsq = af1ch2filt.^2;

z11 = mean(af1ch1filtsq);
z12 = mean(af1ch2filtsq);

loud = -0.691 + 10*log10(z11+z12);

loudness=loud;

end

