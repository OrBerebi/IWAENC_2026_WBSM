classdef InterruptibleAudioPlayer < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        aDW
        is_playing
        fs
        playChannels
        device
        fadeIn
        fadeOut
        bufferSize
        interrupted
        switched
        sigplay
        loop
    end
    
    methods
        
        function obj = InterruptibleAudioPlayer(fs, playChannels, device, fadeIn, fadeOut, bufferSize, loop)
            if isa(device, 'cell')
                device = device{1};
            end
            BitDepth = '24-bit integer';
            obj.fs = fs;
            obj.playChannels = playChannels;
            obj.device = device;
            obj.fadeIn = fadeIn;
            obj.fadeOut = fadeOut;
            obj.bufferSize = bufferSize;    
            obj.loop = loop;
            obj.aDW = audioDeviceWriter('SampleRate',fs,'BitDepth',BitDepth,'Device',device);
%             if ispc
%                 try
%                     obj.aDW.Driver = 'ASIO';
%                 catch
%                 end
%             end
            obj.aDW.ChannelMappingSource = 'Property';
            obj.aDW.ChannelMapping = playChannels;
            obj.is_playing = 0;
            obj.interrupted = 0;
            obj.switched = 0;
            obj.sigplay = [];
        end
        
        function stop(obj)
            if obj.is_playing
                obj.interrupted = 1;
%                 disp('Fading out...');
            else
%                 disp('Not playing!');
            end
        end
        
        function r = isplaying(obj)
            r = obj.is_playing;
%             fprintf('Query: isplaying = %d\n',r);
        end
        
        function changeChannels(obj,playChannels)
            release(obj.aDW);
            obj.aDW.ChannelMapping = playChannels;
%             disp('Changing channels')
        end
        
        function ch = getChannels(obj)
            ch = obj.aDW.ChannelMapping;
        end
        
        function finished = play(obj,sigplay)   
            
            finished = 0;
            
            Y = sigplay;
            L = length(Y);
            
%             d = datetime; d.Format = 'HH:mm:ss.SSS';
%             fprintf('%s: Length is %d samples\n',string(d),L)

            % Apply fade-in
            if obj.fadeIn > 0
                slopelen = obj.fadeIn*obj.fs; % number of samples of the fadein slope
                t = (0:(slopelen-1))';
                slope = 0.5*(1 + cos(pi + 2*pi*t/(2*slopelen))); % raised cosine slope
                win = [slope; ones(L-slopelen,1)];
                win = repmat(win,1,size(Y,2)); % repeat for all channels
                Y = Y.*win;
            end

            % Apply fade-out
            if obj.fadeOut > 0
                slopelen = obj.fadeOut*obj.fs; % number of samples of the fadein slope
                t = (0:(slopelen-1))';
                slope = 0.5*(1 + cos(2*pi*t/(2*slopelen))); % raised cosine slope
                win = [ones(L-slopelen,1);slope];
                win = repmat(win,1,size(Y,2)); % repeat for all channels
                Y = Y.*win;
            end
            
%             d = datetime; d.Format = 'HH:mm:ss.SSS';
%             fprintf('%s: Assigning new signal (rms = %0.8f)\n',string(d),mean(rms(Y)))
            obj.sigplay = Y; % with fade in/out
            
            if ~obj.is_playing
                
%                 disp('Was not playing. Starting to play...');
                
                obj.is_playing = 1;                
                buffer_start = 1;
                stopped = 0;

                while buffer_start <= L

                    drawnow(); % to allow interruption from outside

                    if obj.interrupted
%                         fprintf('Was interrupted. Fading out...\n')
                        obj.interrupted = 0;
                        stopped = 1;
                        remaining_samples = L-buffer_start;
                        slopelen = obj.fadeOut*obj.fs; % number of samples of the fadeout slope
                        if remaining_samples >= slopelen
                            t = (0:(slopelen-1))';
                            slope = 0.5*(1 + cos(2*pi*t/(2*slopelen))); % raised cosine slope
                            slope = repmat(slope,1,size(Y,2)); % repeat for all channels
                            Y(buffer_start:buffer_start+slopelen-1,:) = Y(buffer_start:buffer_start+slopelen-1,:).*slope;
                            Y(buffer_start+slopelen:end,:) = [];
                            L = length(Y);
                        end
                    elseif obj.switched % to switch stimuli
                        obj.switched = 0;
                        % First we apply fade out, just like above
                        remaining_samples = L-buffer_start;
                        slopelen = obj.fadeOut*obj.fs; % number of samples of the fadeout slope
                        if remaining_samples >= slopelen
                            t = (0:(slopelen-1))';
                            slope = 0.5*(1 + cos(2*pi*t/(2*slopelen))); % raised cosine slope
                            slope = repmat(slope,1,size(Y,2)); % repeat for all channels
                            Y(buffer_start:buffer_start+slopelen-1,:) = Y(buffer_start:buffer_start+slopelen-1,:).*slope;
                        end
%                         d = datetime; d.Format = 'HH:mm:ss.SSS';
%                         fprintf('%s: Stimulus was switched. Reminder of signal is swapped with new stimulus (rms = %0.8f)\n',string(d),mean(rms(obj.sigplay)))
                        % Then, we swap the rest of the signal with the new
                        % one
                        ind = buffer_start+slopelen;
                        Y = [Y(1:min(ind-1,end),:);obj.sigplay(min(ind,end):end,:)];
%                         d = datetime; d.Format = 'HH:mm:ss.SSS';
%                         fprintf('%s: New signal is %d samples long. The first %d are from the old signal. The next %d are from the new signal.\n',string(d),length(Y),ind2,length(obj.sigplay)-ind3)
                        L = length(Y);
%                         d = datetime; d.Format = 'HH:mm:ss.SSS';
%                         fprintf('%s: New length is %d samples\n',string(d),L)
                        % Apply fade-in
                        if obj.fadeIn > 0
                            slopelen = obj.fadeIn*obj.fs; % number of samples of the fadein slope
                            t = (0:(slopelen-1))';
                            slope = 0.5*(1 + cos(pi + 2*pi*t/(2*slopelen))); % raised cosine slope
                            slope = repmat(slope,1,size(sigplay,2)); % repeat for all channels
                            ind2 = ind+slopelen;
                            range = ind:min(ind2-1,L);
                            Y(range,:) = Y(range,:).*slope(1:length(range),:);
                        end
                    end

                    audioToDevice = zeros(obj.bufferSize,size(Y,2));

                    buffer_end = buffer_start + obj.bufferSize - 1;

                    if buffer_end >= L
                        if ~obj.loop || stopped
%                             disp('Reached end of stimulus. Stopping...');
                            buffer_end = L;
                            audioToDevice(1:buffer_end-buffer_start+1,:) = Y(buffer_start:buffer_end,:);
                        else
%                             d = datetime; d.Format = 'HH:mm:ss.SSS';
%                             fprintf('%s: Reached end of stimulus and looping. New signal has rms = %0.8f\n',string(d),mean(rms(obj.sigplay)));
                            buffer_end1 = L;
                            buffer_end2 = buffer_end-L;
                            audioToDevice(1:buffer_end1-buffer_start+1,:) = Y(buffer_start:buffer_end1,:);
                            
                            Y = obj.sigplay;
                            L = length(Y);

                            % Apply fade-in. EDIT: already applied to obj.sigplay
%                             if obj.fadeIn > 0
%                                 slopelen = obj.fadeIn*obj.fs; % number of samples of the fadein slope
%                                 t = (0:(slopelen-1))';
%                                 slope = 0.5*(1 + cos(pi + 2*pi*t/(2*slopelen))); % raised cosine slope
%                                 win = [slope; ones(L-slopelen,1)];
%                                 win = repmat(win,1,size(Y,2)); % repeat for all channels
%                                 Y = Y.*win;
%                             end
                            
                            audioToDevice(buffer_end1-buffer_start+2:end,:) = Y(1:buffer_end2,:);
                        end
                    else
%                         fprintf('%d - %d (L = %d)\n',buffer_start,buffer_end,L)
                        audioToDevice = Y(buffer_start:buffer_end,:);
                    end

                    numUnderrun = obj.aDW(audioToDevice);

                    if numUnderrun > 0
%                         fprintf('Between t=%f and t=%f, Audio player queue was underrun by %d samples.\n',buffer_start/obj.fs,buffer_end/obj.fs,numUnderrun);
                        sync = 0;
                    end

                    buffer_start = buffer_start + obj.bufferSize;
                    if obj.loop && buffer_start > L
                        if ~stopped
                            buffer_start = buffer_start-L;
                        else
                            stopped = 0;
                        end
                    end
                        
                end
%                 disp('Finished playing')
                obj.is_playing = 0;
                finished = 1;
                
            else
%                 d = datetime; d.Format = 'HH:mm:ss.SSS';
%                 fprintf('%s: Already playing. Activating Switched flag...\n',string(d))              
                obj.switched = 1;
                
            end

        end
        
        function delete(obj)
            release(obj.aDW);
        end
        
    end
end

