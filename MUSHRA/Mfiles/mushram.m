function mushram(varargin)

% MUSHRAM   Runs the MUSHRA experiments described in the configuration file
% mushram_config.txt
%
% mushram runs both training and evaluation phases
%
% mushram('training') runs the training phase, followed by training2 and
% evaluation
%
% mushram('training2') runs the second training phase (familiarization),
% followed by evaluation
%
% mushram('evaluation') runs the evaluation phase only
%
% mushram(phase, initPos) where initPos is the initial position of the GUI
% window. Can be used to make sure the GUI appears where is wanted
%
% mushram(phase, [], 'no_random') not randomize the order of the experiments (but does randomize the
% order of the test files within each experiment anyway)
%
% mushram(phase, [], expe_order) performs the evaluation experiments in the
% alternative order given by the vector expe_order (must span all the
% experiments defined in the configuration file)
%
% mushram(phase, [], [], randomFiles) - default: true. if false, not
% randomize the order of the test files
%
% Based on: http://c4dm.eecs.qmul.ac.uk/downloads/#mushram
% Edited by: Zamir Ben-Hur
% Last updated: July 2021
%

%%%process input flags
global nowPlaying
global mushra_settings
global played
global IR_FLAG


%default parameters
phase=0;
random_expe=true;
initPos = [];

%user parameters
myPath = fileparts(mfilename('fullpath'));
if ~isempty(varargin)
    phasename=varargin{1};
    
    if ischar(phasename)
        if strcmp(phasename,'training')
            random_expe = true;
            if IR_FLAG
                configFile = fullfile(myPath,'..','ConfigFiles','MUSHRA_config_training_IRs.txt');
            else
                configFile = fullfile(myPath,'..','ConfigFiles','MUSHRA_config_training.txt');
            end
        elseif strcmp(phasename,'training2')
            phase=2;
            randomFiles = true;
            if IR_FLAG
                configFile = fullfile(myPath,'..','ConfigFiles','MUSHRA_config_training2_IRs.txt');
            else
                configFile = fullfile(myPath,'..','ConfigFiles','MUSHRA_config_training2.txt');
            end
        elseif strcmp(phasename,'evaluation')
            phase=1;
            if IR_FLAG
                configFile = fullfile(myPath,'..','ConfigFiles','MUSHRA_config_IRs.txt');
            else
                configFile = fullfile(myPath,'..','ConfigFiles','MUSHRA_config.txt');
            end
        elseif ~strcmp(phasename,'all')
            errordlg('Bad input parameters.','Error');
            return;
        end
    else
        errordlg('Bad input parameters.','Error');
        return;
    end
    if length(varargin) > 1
        initPos = varargin{2};
    end
    
    if length(varargin) > 2
        norandom=varargin{3};
        if ischar(norandom)
            if strcmp(varargin{4},'no_random')
                random_expe=0;
            else
                errordlg('Bad input parameters.','Error');
                return;
            end
        elseif isempty(norandom)
            
        elseif isnumeric(norandom)
            expe_order=norandom;
        else
            errordlg('Bad input parameters.','Error');
            return;
        end
    end
    if length(varargin) > 3
        randomFiles=varargin{4};
    end
end

nowPlaying = [];
mushra_settings.done = 0;

% Target loudness (LUFS)
% Calibrated for RME settings: headphone gain = -5dB, speaker gain = -5dB
targetLoudness = -6; %-36;

%% parsing the config file
%opening the file
fid=fopen(configFile,'r');
if fid~=-1
    config=fscanf(fid,'%c');
    fclose(fid);
else
    errordlg(['The configuration file ' configFile ' does not exist or cannot be opened.'],'Error');
    return;
end

%suppressing spurious linefeeds and transforming Windows linefeeds into Unix ones
if length(config)>1
    config=strrep(config,char(13),'');
    c=find(config~=10);
    config=[config(1:c(end)) newline newline];
    while contains(config,[newline newline newline])
        config=strrep(config,[newline newline newline],[newline newline]);
    end
else
    errordlg(['The configuration file ' configFile ' is empty.'],'Error');
    return;
end


%parsing and checking the file names for all experiments
dblines=strfind(config,[newline newline]);
nbexpe=length(dblines); % handles.NumStim
% expconfig=config(1:dblines(1));
% lines=strfind(expconfig,newline);
nbfile=zeros(nbexpe,1);
% files=cell(nbexpe,nbfile);
stimFiles = cell(nbexpe,1);
dblines=[-1 dblines];
if IR_FLAG
    preLines = 2;
else
    preLines = 1;
end
for e = 1:length(dblines)-1
    expconfig=config(dblines(e)+2:dblines(e+1));
    lines=strfind(expconfig,newline);
    nbfile(e) = length(lines)-preLines;
    %     if length(lines) == nbfile+1
    lines=[0 lines];
    for f=1:length(lines)-1
        if f==1 % metric
            metricString{e} = expconfig(lines(f)+1:lines(f+1)-1);
        elseif IR_FLAG && f==2 % folder of stimulus
            stimFile = expconfig(lines(f)+1:lines(f+1)-1);
            if exist(fullfile(myPath,'..',stimFile),'file')
                stimFiles{e}=stimFile;
            else
                errordlg(['The specified sound file ' stimFile ' does not exist. Check the configuration file ' configFile '.'],'Error');
                return;
            end
        else
            
            file = expconfig(lines(f)+1:lines(f+1)-1);
            if exist(fullfile(myPath,'..',file),'file')
                files{e,f-preLines}=file;
            else
                errordlg(['The specified sound file ' file ' does not exist. Check the configuration file ' configFile '.'],'Error');
                return;
            end
            [~,name,~] = fileparts(file);
            
            name = regexprep(name, ' ', '_');
            
            index_underscore = regexp(name,'_');
            if length(index_underscore)>4
                mushra_settings.responses.stimuli{e,f-preLines} = name([1:index_underscore(1) index_underscore(3)+1:index_underscore(4) index_underscore(5)+1:end ]); %regexprep(name(1:index_underscore(1)-1), ' ', '_');
            else
                mushra_settings.responses.stimuli{e,f-preLines} = name;
                
            end
            if IR_FLAG
                mushra_settings.responses.stimuli{e,f-preLines} = regexprep(mushra_settings.responses.stimuli{e,f-preLines},'-','');
                [~,stimFolderSplit,~] = fileparts(stimFiles{e});
                mushra_settings.responses.stimuli{e,f-preLines} = [ mushra_settings.responses.stimuli{e,f-preLines} '_' stimFolderSplit];
            end
        end
    end
    %     else
    %         errordlg(['The number of test files must be the same for all experiments. Check the configuration file ' configFile '.'],'Error');
    %         return;
    %     end
    
end
mushra_settings.stimFiles = stimFiles;
%%

%%%randomizing the order of the experiments or checking the alternative experiment order is acceptable
if exist('expe_order')
    err=false;
    for e=1:nbexpe
        if ~any(expe_order==e)
            err=true;
        end
    end
    if err
        errordlg('Bad input parameters.','Error');
        return;
    end
elseif random_expe
    expe_order=randperm(nbexpe);
else
    expe_order=1:nbexpe;
end

switch phase
    case 1
        %%%opening the GUI for the evaluation phase
        %asking for the name of the results file
        dlisting = dir(fullfile(myPath,'..','Results','*.csv'));
        id = length(dlisting)+1;
        mushra_settings.responses.id = num2str(id,'%02d');
        if strcmp(mushra_settings.hpcf.hpName,'Other')
            mushra_settings.responses.id = [mushra_settings.responses.id '_noHpEq'];
        end
        mushra_settings.responses.date = datestr(clock);
        %     [filename,pathname] = uiputfile('mushram_results.txt','Results file name');
        %     resultfile = [pathname filename];
        %     if ~resultfile
        %         return;
        %     end
        mushra_settings.resultfile = fullfile(myPath,'..','Results',['MUSHRA_Results_' mushra_settings.responses.userData.name mushra_settings.responses.id '.csv']);
        if ~exist(fullfile(myPath,'..','Results'),'dir')
            mkdir(fullfile(myPath,'..','Results'))
        end
        %opening the GUI
        
        fig = evaluation_gui(nbfile(expe_order(1)), 1, nbexpe, initPos, metricString{expe_order(1)});
        if isempty(fig)
            errordlg('There are too many test files to display. Try increasing the resolution of the screen.','Error');
            return;
        end
        
        handles=guihandles(fig);
        handles.fig = fig;
        
        loadFiles(stimFiles, nbexpe, nbfile, files , targetLoudness);
        
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        userPrompt = sprintf(['\n                                       Elvaluation Phase:\n\n- Rate the similarity between the test signals and the reference with respect to ' metricString{expe_order(1)} '.\n- You must give 100 points to the hidden reference, and rank the other test signals in order that reflects the audible impairments.\n- You can rate seneral signals with 100.\nIMPORTANT: If you do not hear a diffrence between several signals, please give them that same rating.']);
        h=msgbox(['\fontsize{16}' userPrompt],'Elvaluation',CreateStruct);
        uiwait(h);
        
        %randomizing the order of the tested files
        file_order=randperm(nbfile(expe_order(1)));
        
        
    case 0
        randomFiles = true;
        %%%opening the GUI for the training phase
        %opening the GUI
        fig = training_gui(nbfile(expe_order(1)), 1, nbexpe, initPos, metricString{expe_order(1)});
        if isempty(fig)
            errordlg('There are too many experiments or too many test files to display. Try increasing the resolution of the screen.','Error');
            return;
        end
        handles=guihandles(fig);
        handles.fig = fig;
        
        loadFiles(stimFiles, nbexpe, nbfile, files , targetLoudness);
        
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        userPrompt = sprintf(['\n                                       User Interface Training:\n\n- Rate the similarity between the test signals and the reference with respect to ' metricString{expe_order(1)} '.\n- You must give 100 points to the hidden reference, and rank the other test signals in order that reflects the audible impairments.\n']);

        h=msgbox(['\fontsize{16}' userPrompt],'Training',CreateStruct);
        uiwait(h);
        
        if randomFiles
            %randomizing the order of the tested files
            file_order=randperm(nbfile(expe_order(1)));
        else
            file_order = repmat(1:nbfile(expe_order(1)),nbexpe,1);
        end
        
        
        rgbImage = imread(fullfile(myPath,'successIcon.png'));
        mushra_settings.myIcon = imresize(rgbImage, [64, 64]);
        
    case 2 % familiarization training
        mushra_settings.randomFiles = randomFiles;
        fig = training2_gui(nbfile, nbexpe, 1, initPos, metricString);
        handles=guihandles(fig);
        handles.fig = fig;
        loadFiles(stimFiles, nbexpe, nbfile, files , targetLoudness);
        
        if randomFiles
            CreateStruct.Interpreter = 'tex';
            CreateStruct.WindowStyle = 'modal';
            userPrompt = sprintf('\n                                   Signals Familiarization:\n\n- These are the test signals that you will evaluate next.\n- In this stage you will be presented with all test signals that you will evaluate in the listening test.\n- The goal of this stage is to expose you to the full range and nature of impairments that you will experience during the test.\n- You may take some time now to become familiar with the signals.\n- Note that the signals are listed in a random order.\n- You must listen at least once to each signal.\n');
            h=msgbox(['\fontsize{16}' userPrompt],'Training',CreateStruct);
            uiwait(h);
        end
        for e=1:nbexpe
            if randomFiles
                %randomizing the order of the tested files
                file_order_tmp = randperm(nbfile(e));
                file_order_tmp(file_order_tmp==1) = [];
                file_order(e,:) = [1 file_order_tmp]; % in order to alway have the reference first
            else
                file_order(e,:) = 1:nbfile(e);
            end
            
        end
end

played = zeros(nbexpe,max(nbfile));

devinfo = getAudioDevInfo_onlyWriter;
% mushra_settings.SoundcardID = 1;
mushra_settings.device = devinfo(mushra_settings.SoundcardID,1);
fadeIn  = 0.001; % fade-in [sec]
fadeOut = 0.001; % fade-out [sec]
bufferSize = 1024; % buffer size [samples] (for example, 2048; 1024; 256; 512)
loop = true;
% handles.headphoneCh = [1 2];
mushra_settings.player = InterruptibleAudioPlayer(mushra_settings.fs,mushra_settings.headphoneCh,mushra_settings.device, fadeIn, fadeOut, bufferSize, loop);


%%%storing data within the GUI
mushra_settings.expe=1;
mushra_settings.ratings=zeros(nbexpe,max(nbfile));
mushra_settings.expe_order=expe_order;
mushra_settings.nbexpe=nbexpe;
mushra_settings.file_order=file_order;
mushra_settings.nbfile=nbfile;
mushra_settings.files=files;
mushra_settings.time=clock;
mushra_settings.nowPlaying = [];
mushra_settings.metricString = metricString;

guidata(handles.fig,handles);


function loadFiles(stimFiles, nbexpe, nbfile, files, targetLoudness )
% load files
global mushra_settings
global IR_FLAG
myPath = fileparts(mfilename('fullpath'));
fbar = waitbar(0,'Loading, please wait...','WindowStyle','modal');
k=1;
for e=1:nbexpe
    if IR_FLAG
        try
            [yStim,fs]=audioread(fullfile(myPath,'..',stimFiles{e}));
        catch ME
            rethrow(ME)
        end
        
        if e==1
            mushra_settings.fs = fs;
            if ~strcmp(mushra_settings.hpcf.hpName,'Other')
                if mushra_settings.hpcf.fs ~= mushra_settings.fs
                    warning('Headphone filter and data sample rate mismatch!. Resample HpEq.')
                    %                     mushra_settings.hpcf = [];
                    %                     mushra_settings.hpcf.hpName = 'Other';
                    [p,q]=rat(mushra_settings.fs/ mushra_settings.hpcf.fs);
                    mushra_settings.hpcf.minPhase = resample(mushra_settings.hpcf.minPhase,p,q);
                end
            end
        elseif fs~=mushra_settings.fs
            error(['Sample frequency of file ' stimFiles{e} ' is not ' num2str(mushra_settings.fs)])
        end
    end
    
    loudness = zeros(length(nbfile(e)),1);
    stim_len = zeros(length(nbfile(e)),1);
    for f=1:nbfile(e)
        waitbar(k/(nbexpe*max(nbfile)),fbar);
        [y,fs]=audioread(fullfile(myPath,'..',files{e,f}));
        if IR_FLAG
            if fs~=mushra_settings.fs
                [p,q]=rat(mushra_settings.fs/fs);
                y = resample(y,p,q);
                
                %             error(['Sample frequency of file ' files{e,f} ' is not ' num2str(mushra_settings.fs)])
            end
        else
            if e==1 && f==1
                mushra_settings.fs = fs;
                if ~strcmp(mushra_settings.hpcf.hpName,'Other')
                    if mushra_settings.hpcf.fs ~= mushra_settings.fs
                        warning('Headphone filter and data sample rate mismatch!. Resample HpEq.')
                        %                     mushra_settings.hpcf = [];
                        %                     mushra_settings.hpcf.hpName = 'Other';
                        [p,q]=rat(mushra_settings.fs/ mushra_settings.hpcf.fs);
                        mushra_settings.hpcf.minPhase = resample(mushra_settings.hpcf.minPhase,p,q);
                    end
                end
            else
                if fs~=mushra_settings.fs
                    [p,q]=rat(mushra_settings.fs/fs);
                    y = resample(y,p,q);
                    
                    %             error(['Sample frequency of file ' files{e,f} ' is not ' num2str(mushra_settings.fs)])
                end
            end
        end
        
        if ~strcmp(mushra_settings.hpcf.hpName,'Other')
            y = miro_fastConv(y,mushra_settings.hpcf.minPhase);
        end
        
        % if IR_FLAG is true, convolve IR with stimulus
        if IR_FLAG
            y = miro_fastConv(y,yStim);
        end
        
        % if mono - transform to stereo
        if size(y,2)==1
            y = [y y];
        end
        % remove zeros
        %         [row_1,~] = find(abs(y)>1e-10,1);
        %         [row_last,~] = find(abs(y)>1e-10,1,'last');
        %         y = y(row_1:row_last,:);
        
        loudness(f) = integratedLoudness(y,mushra_settings.fs);
        stim_len(f) = size(y,1);
        
        mushra_settings.players{e,f} = y; %audioplayer(y,handles.fs);
        %         handles.players{e,f}.StopFcn = {@player_stop};
        k=k+1;
    end
    
    % Zeropad if needed
    maxlen = max(stim_len);
    for ind_cond = 1:nbfile(e)
        delta = maxlen - stim_len(ind_cond);
        mushra_settings.players{e,ind_cond} = [mushra_settings.players{e,ind_cond};zeros(delta,2)];
    end
    
    % Normalize perceived loudness to some target value (LUFS)
    for ind_cond = 1:nbfile(e)
        delta = db2mag(targetLoudness - loudness(ind_cond));
        mushra_settings.players{e,ind_cond} = mushra_settings.players{e,ind_cond}*delta;
    end
    
    % Prevent clipping and warn if it occurs
    sigmax_all = 0;
    for ind_cond = 1:nbfile(e)
        sigmax = max(max(abs(mushra_settings.players{e,ind_cond})));
        %            if sigmax > 1
        %                warning('Signal clipped in condition %s, program "%s" (sigmax = %f)',cond.Name,prog.Name,sigmax)
        %            end
        if sigmax > sigmax_all
            sigmax_all = sigmax;
        end
    end
    %    fprintf('Peak for exp %d = %0.2f dBFS\n',e,db(sigmax_all));
    if sigmax_all > 1
        %        warning("Clipping occured! Normalizing exp %d by %0.2f dB...",e,db(sigmax_all))
        for ind_cond = 1:nbfile(e)
            mushra_settings.players{e,ind_cond} = mushra_settings.players{e,ind_cond}./sigmax_all;
        end
    end
    
end
close(fbar)


function ab = miro_fastConv(a,b)

% Internal use only

NFFT = size(a,1)+size(b,1)-1;
A    = fft(a,NFFT);
B    = fft(b,NFFT);
AB   = A.*B;
ab   = ifft(AB);


function audioData2 = loudnessEq(audioData1,audioData2)

threshold = 1;
errorThreshold = threshold / 100;


% Compute loudness
loudness_ref = LoudnessCalc_ITU(audioData1);
loudness_2 = LoudnessCalc_ITU(audioData2);

% start with stepsize of 10 dB
if loudness_ref > loudness_2
    ampStep = 10;
else
    ampStep = -10;
end

lastStatus = loudness_ref > loudness_2;
% break if too many iterations are needed
whileCounter = 0;
whileLimit = 25;

currentAmp = 0;
multFactor = 1;

while (abs(abs(loudness_ref-loudness_2) / loudness_ref) > errorThreshold)
    
    if lastStatus ~= ( loudness_ref > loudness_2)    % if status change => decrease and invert stepsize
        ampStep = ampStep / -10;
        lastStatus = loudness_ref > loudness_2;
    end
    
    currentAmp = currentAmp + ampStep;
    
    multFactor = 10.^(currentAmp./20);
    
    loudness_2 = LoudnessCalc_ITU(multFactor.*audioData2);
    
    
    whileCounter = whileCounter +1;
    if (whileCounter >= whileLimit)
        warning('Difference threshold can not be met. Try increasing it');
        return
    end
end
audioData2 = multFactor.*audioData2;