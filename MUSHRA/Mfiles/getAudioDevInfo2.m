function audioInfo = getAudioDevInfo2()
%GETAUDIODEVINFO2 returns a llist of recognized full-duplex audio devices.
% As a general rule, it will recognize ASIO (Windows) and CoreAudio
% (macOS), but not generic Windows drivers such as DirectSound.
%
% *requires "Audio System Toolbox"
% 
%USAGE
%   audioInfo = getAudioDevInfo2;
%
%INPUT PARAMETERS
%   none
%
%OUTPUT PARAMETERS
%   audioInfo : cell with all full-duplex audio devices and their
%               respective maximum number of input/ouput channels
% 
%
%Author: Isaac Engel (August, 2018)

audioInfo = cell(1,2);

aPR = audioDeviceWriter;

dev_list = getAudioDevices(aPR);
num_in_dev = numel(dev_list);

for i = 1:num_in_dev
    
    audioInfo{i,1} = [dev_list{i}];

    aPR = audioDeviceWriter('Device',dev_list{i});
    
    tmp_dev_info = info(aPR);
    audioInfo{i,2} = tmp_dev_info.MaximumRecorderChannels;
    audioInfo{i,3} = tmp_dev_info.MaximumPlayerChannels;
end

release(aPR);
delete(aPR);
