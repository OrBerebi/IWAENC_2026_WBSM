%%  
% MUSHRA main script
% 


clear all
close all force
clc
set(0, 'DefaultFigureWindowStyle', 'normal');

restoredefaultpath

% Add Main Path
myPath = fileparts(mfilename('fullpath')); 
addpath(fullfile(myPath, 'Mfiles'));

global IR_FLAG


%%

skipIntro = false; % flag to skip introduction screen and info collection screen (only choose Headphone type) 

IR_FLAG = false; % if true - using IRs + dry audio with convolution during the test, else, use wav files as they are


% Run exp
if skipIntro
    MUSHRA_training_App_noIntro
else
    if verLessThan('matlab','9.6')
        MUSHRA_training_App_2017a
    else
        MUSHRA_training_App
    end
end