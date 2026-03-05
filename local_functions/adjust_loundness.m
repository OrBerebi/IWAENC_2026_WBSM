function y = adjust_loundness(x,fs,targetLUFS)
% Compute LUFS loudness of each audio signal
%channelWeights = ones([1 size(x,2)]);
%loudness1 = integratedLoudness(x, fs,   channelWeights);
loudness1 = integratedLoudness(x, fs);
% Calculate the gain adjustment in dB to match targetLUFS
gain1 = targetLUFS - loudness1;
% Apply gain adjustment
y = x * 10^(gain1 / 20);
end