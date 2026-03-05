function audioInfo = getAudioDevInfo_onlyWriter()
%GETAUDIODEVINFO2 returns a llist of recognized writer audio devices.
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
%   audioInfo : cell with all audio devices and their
%               respective maximum number of ouput channels
% 
%
%Author: Isaac Engel (August, 2018)
%Edited by Zamir Ben-Hur (May 2020)

audioInfo = cell(1,2);

aPR = audioDeviceWriter;

dev_list = getAudioDevices(aPR);
num_in_dev = numel(dev_list);

for i = 1:num_in_dev
    
    audioInfo{i,1} = [dev_list{i}];

    aPR = audioDeviceWriter('Device',dev_list{i});
    
    tmp_dev_info = info(aPR);
    audioInfo{i,2} = tmp_dev_info.MaximumOutputChannels;
end

release(aPR);
delete(aPR);
