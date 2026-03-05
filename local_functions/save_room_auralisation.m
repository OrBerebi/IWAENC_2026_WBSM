function save_room_auralisation(res_out,fs,save_path,sig_path,brir_flag)

targetLUFS = -35;



[~, name, ext] = fileparts(sig_path);
stim_name = name;  % Combines name and extension
save_dir = save_path +"/room_auralisation/" + stim_name;
mkdir(save_dir)


[s, signal_fs] = audioread(sig_path);
if (signal_fs ~= fs)
    s = resample(s, fs, signal_fs);
    signal_fs = fs;
end
 
filename = save_dir + "/ref.wav";
BRIR     = res_out.BRIR_ref;
if ~brir_flag
    y   = fftfilt(BRIR, s);
    y   = adjust_loundness(y,fs,targetLUFS);
else
    y   = BRIR;
    y   = y./max(max(abs(y)));
    filename = save_dir + "/BRIR ref.wav";
end
audiowrite(filename,y,fs);


filename = save_dir + "/FOV-BSM-ls.wav";
BRIR     = res_out.BRIR_ls_fov;
if ~brir_flag
    y   = fftfilt(BRIR, s);
    y   = adjust_loundness(y,fs,targetLUFS);
else
    y   = BRIR;
    y   = y./max(max(abs(y)));
    filename = save_dir + "/BRIR FOV-BSM-ls.wav";
end
audiowrite(filename,y,fs);

filename = save_dir + "/FOV-BSM-mls.wav";
BRIR     = res_out.BRIR_mls_fov;
if ~brir_flag
    y   = fftfilt(BRIR, s);
    y   = adjust_loundness(y,fs,targetLUFS);
else
    y   = BRIR;
    y   = y./max(max(abs(y)));
    filename = save_dir + "/BRIR FOV-BSM-mls.wav";
end
audiowrite(filename,y,fs);

filename = save_dir + "/BSM-ls.wav";
BRIR     = res_out.BRIR_ls_normal;
if ~brir_flag
    y   = fftfilt(BRIR, s);
    y   = adjust_loundness(y,fs,targetLUFS);
else
    y   = BRIR;
    y   = y./max(max(abs(y)));
    filename = save_dir + "/BRIR BSM-ls.wav";
end
audiowrite(filename,y,fs);

filename = save_dir + "/BSM-mls.wav";
BRIR     = res_out.BRIR_mls_normal;
if ~brir_flag
    y   = fftfilt(BRIR, s);
    y   = adjust_loundness(y,fs,targetLUFS);
else
    y   = BRIR;
    y   = y./max(max(abs(y)));
    filename = save_dir + "/BRIR BSM-mls.wav";
end
audiowrite(filename,y,fs);



end

function y = adjust_loundness(x,fs,targetLUFS)

% Compute LUFS loudness of each audio signal
loudness1 = integratedLoudness(x, fs);
% Calculate the gain adjustment in dB to match targetLUFS
gain1 = targetLUFS - loudness1;
% Apply gain adjustment
y = x * 10^(gain1 / 20);
end