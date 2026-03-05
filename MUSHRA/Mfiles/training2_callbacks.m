function training2_callbacks(varargin)

% training_CALLBACKS   Callback functions for the training interface
%
% training_callbacks(fname,varargin) executes the callback function
% fname_callback with various parameters

fname=[varargin{1},'_callback'];
feval(fname,varargin{2:end});



%%%saving the rating results and proceeding to the next experiment or exiting
function proceed_callback(handles)
global nowPlaying played mushra_settings
totalStim = 0;
for e=1:mushra_settings.nbexpe
    totalStim =  totalStim + mushra_settings.nbfile(e);
end
if ~(sum(sum(played)) == totalStim)
    CreateStruct.Interpreter = 'tex';
    CreateStruct.WindowStyle = 'modal';
    userPrompt = sprintf('\nYou have to listen to all sounds before proceeding.\n');
    h=msgbox(['\fontsize{16}' userPrompt], 'Warning', 'warn',CreateStruct);
    uiwait(h);
    return
end
if ~isempty(nowPlaying)
    unhighlight_button(handles,nowPlaying);
end
nowPlaying = [];
stop_callback(handles)

% done training, continue to evaluation
initPos = handles.fig.Position(1:2);
mushra_settings.done = 1;
close(handles.fig);
mushram('evaluation',initPos);



%%%stop sound files
function stop_callback(handles)
global nowPlaying mushra_settings
stop(mushra_settings.player)
if ~isempty(nowPlaying)
    unhighlight_button(handles,nowPlaying);
end
nowPlaying = [];


%%%playing sound files
function play_callback(handles,e,f)
global nowPlaying mushra_settings played
if ~isempty(nowPlaying)
    unhighlight_button(handles,nowPlaying);
end
if f
    playNum = mushra_settings.file_order(e,f+1);
    nowPlaying = [e f];
    played(e,f+1) = 1;
else
    playNum = 1;
    nowPlaying = [e 0];
    played(e,1) = 1;
end


highlight_button(handles,nowPlaying);


finishedPlaying = play(mushra_settings.player,mushra_settings.players{e,playNum});
if finishedPlaying && ~isempty(nowPlaying)
    unhighlight_button(handles,nowPlaying);
    nowPlaying = [];
end

function highlight_button(handles,button)
eval(sprintf("set(handles.play%d_%d,'backgroundcolor',[0.3,1,0.3]);",button(1),button(2)));


function unhighlight_button(handles,button)
global played
if played(button(1),button(2)+1)
    eval(sprintf("set(handles.play%d_%d,'backgroundcolor',[0.6,0.8,1]);",button(1),button(2)));
else
    eval(sprintf("set(handles.play%d_%d,'backgroundcolor',[0.94 0.94 0.94]);",button(1),button(2)));
end


function keyPressCallback_callback(hObject, eventdata, handles)
global mushra_settings nowPlaying
switch eventdata.Key
    case 'q'
        stop_callback(handles)
    case 'return'
        proceed_callback(handles)
    case 'w'
        if  ~isempty(nowPlaying) && nowPlaying(1) > 1
            play_callback(handles,nowPlaying(1)-1, nowPlaying(2))
        end
    case 's'
        if  ~isempty(nowPlaying) && nowPlaying(1) < mushra_settings.nbexpe
            play_callback(handles,nowPlaying(1)+1,nowPlaying(2))
        end
    case 'a'
        if  ~isempty(nowPlaying) && nowPlaying(2)>0
            play_callback(handles,nowPlaying(1),nowPlaying(2)-1)
        end
        
    case 'd'
        if isempty(nowPlaying)
            play_callback(handles,1,0)
        else
            if nowPlaying(2) < mushra_settings.nbfile(nowPlaying(1))-1
                play_callback(handles,nowPlaying(1),nowPlaying(2)+1)
            end
        end
end

function CloseFigure_callback(hObject, eventdata, handles)
stop_callback(handles)
% Close request function
% to display a question dialog box
global mushra_settings
if mushra_settings.done
    delete(gcf)
else
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'modal';
CreateStruct.Default = 'No';
selection = questdlg('\fontsize{16}Are you sure you want to exit the experiment?',...
    'Close Request Function',...
    'Yes','No',CreateStruct);
switch selection
    case 'Yes'
        delete(gcf)
    case 'No'
        return
end
end
