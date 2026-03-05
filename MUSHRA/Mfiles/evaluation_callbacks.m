function evaluation_callbacks(varargin)

% EVALUATION_CALLBACKS   Callback functions for the evaluation interface
%
% evaluation_callbacks(fname,varargin) executes the callback function
% fname_callback with various parameters

fname=[varargin{1},'_callback'];
feval(fname,varargin{2:end});



%%%saving the rating results and proceeding to the next experiment or exiting
function results_callback(handles)
global nowPlaying mushra_settings played
if ~isempty(nowPlaying)
    unhighlight_button(handles,nowPlaying);
    unhighlight_slider(handles,nowPlaying);
end
nowPlaying = [];
stop_callback(handles)


totalStim = mushra_settings.nbfile(mushra_settings.expe_order(mushra_settings.expe))+1;

if ~(sum(played(mushra_settings.expe_order(mushra_settings.expe),:)) == totalStim)
    CreateStruct.Interpreter = 'tex';
    CreateStruct.WindowStyle = 'modal';
    userPrompt = sprintf('\nYou have to listen to all sounds before proceeding.\n');
    h=msgbox(['\fontsize{16}' userPrompt], 'Warning', 'warn',CreateStruct);
    uiwait(h);
    return
end

%getting the ratings for all files
for f=1:mushra_settings.nbfile(mushra_settings.expe_order(mushra_settings.expe))
    mushra_settings.ratings(mushra_settings.expe_order(mushra_settings.expe),mushra_settings.file_order(f))=get(getfield(handles,['slider' int2str(f)]),'Value');
end
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'modal';
if ~any(mushra_settings.ratings(mushra_settings.expe_order(mushra_settings.expe),:)==100)
    userPrompt = sprintf('\nAt least one of the sounds must be rated 100.\n');
    h=msgbox(['\fontsize{16}' userPrompt], 'Warning', 'warn',CreateStruct);
    uiwait(h);
    return;
end

if mushra_settings.expe<mushra_settings.nbexpe
    %saving the whole results (to avoid losing data if the program terminates early)
    saveResults()
    
    mushra_settings.expe=mushra_settings.expe+1;
    if mushra_settings.nbfile(mushra_settings.expe_order(mushra_settings.expe)) ~= mushra_settings.nbfile(mushra_settings.expe_order(mushra_settings.expe-1))
        
        initPos = handles.fig.Position(1:2);
        mushra_settings.done = 1;
        close(handles.fig);
        fig=evaluation_gui(mushra_settings.nbfile(mushra_settings.expe_order(mushra_settings.expe)),mushra_settings.expe,mushra_settings.nbexpe,initPos, mushra_settings.metricString{mushra_settings.expe_order(mushra_settings.expe)});
        handles=guihandles(fig);
        handles.fig = fig;
        nowPlaying = [];
        mushra_settings.done = 0;
    end
    
    if ~strcmp(mushra_settings.metricString{mushra_settings.expe_order(mushra_settings.expe)},mushra_settings.metricString{mushra_settings.expe_order(mushra_settings.expe-1)}) % if not same metric, popup note
        CreateStruct.Interpreter = 'tex';
        CreateStruct.WindowStyle = 'modal';
        userPrompt = sprintf(['\n  Note that the evalution metric is changed to: ' mushra_settings.metricString{mushra_settings.expe_order(mushra_settings.expe)} '.\n']);
        h=warndlg(['\fontsize{16}' userPrompt],'Elvaluation',CreateStruct);
        uiwait(h);
        set(handles.experiment2,'String',['Rate the similarity between the test signals and the reference with respect to ' mushra_settings.metricString{mushra_settings.expe_order(mushra_settings.expe)} '.']);
    end
    
    %updating title
    set(handles.experiment,'String',['Experiment ' int2str(mushra_settings.expe) '/' int2str(mushra_settings.nbexpe)]);
    if mushra_settings.expe==mushra_settings.nbexpe
        pos=get(handles.results,'Position');
        pos(1)=pos(1)+2.5;
        pos(3)=19;
        set(handles.results,'Position',pos,'String','Save & exit (enter)');
    end
    %moving all the sliders back to 0
    for f=1:mushra_settings.nbfile(mushra_settings.expe_order(mushra_settings.expe))
        shandle=getfield(handles,['slider' int2str(f)]);
        set(shandle,'Value',0);
        rhandle=getfield(handles,['rating' int2str(f)]);
        set(rhandle,'String',0);
    end
    %randomizing the order of the tested files for the next experiment
    mushra_settings.file_order=randperm(mushra_settings.nbfile(mushra_settings.expe_order(mushra_settings.expe)));
    %testing whether a break is needed before the next experiment
    if etime(clock,mushra_settings.time) > 20*60
        wfig=warndlg(['You have been working for ' int2str(round(etime(clock,mushra_settings.time)/60)) 'minutes. It is recommended that you take a break of at least the same duration before starting the next experiment. Click on OK when you are ready.'],'Warning');
        uiwait(wfig);
    end
    mushra_settings.time=clock;
    
    
    guidata(handles.fig,handles);
else
    %exiting
    mushra_settings.responses.Comment = commentDialbox;
    
    saveResults()
    mushra_settings.done = 1;
    close(handles.fig);
end




function saveResults()
global mushra_settings
id = [ '_' mushra_settings.responses.userData.name mushra_settings.responses.id];
responses = mushra_settings.responses;
responses.ratings = mushra_settings.ratings;
responses.headphone = mushra_settings.hpcf.hpName;
save([mushra_settings.resultfile(1:end-4) '.mat'],'responses');


out{1,1} = 'id';
out{1,2} = mushra_settings.responses.id;
out{2,1} = 'date';
out{2,2} = mushra_settings.responses.date;
k=3;
for e=1:mushra_settings.nbexpe
    for f=1:mushra_settings.nbfile(e)
        out{k,1} = [mushra_settings.responses.stimuli{e,f},];
        out{k,2} = mushra_settings.ratings(e,f);
        k=k+1;
    end
end

% Save to csv file
% Convert cell to a table and use first row as variable names
T = cell2table(out);

% Write the table to a CSV file
writetable(T,mushra_settings.resultfile)


%%%rounding and displaying the values of the sliders
function slider_callback(handles,f)

shandle=getfield(handles,['slider' int2str(f)]);
set(shandle,'Value',round(get(shandle,'Value')));
rhandle=getfield(handles,['rating' int2str(f)]);
set(rhandle,'String',get(shandle,'Value'));


%%%stop sound files
function stop_callback(handles)
global nowPlaying mushra_settings
% for e=1:mushra_settings.nbexpe
%     for f=1:mushra_settings.nbfile(e)
%         mushra_settings.players{e,f}.StopFcn = [];
%         stop(mushra_settings.players{e,f});
%     end
% end
stop(mushra_settings.player)
if ~isempty(nowPlaying)
    unhighlight_button(handles,nowPlaying);
    unhighlight_slider(handles,nowPlaying);
end
nowPlaying = [];

% guidata(handles.fig,handles);

%%%playing sound files
function play_callback(handles,f)
global nowPlaying mushra_settings played
% currentSamp = 1;
if ~isempty(nowPlaying)
    unhighlight_button(handles,nowPlaying);
    unhighlight_slider(handles,nowPlaying);
end
if f
    playNum = mushra_settings.file_order(f);
    nowPlaying = f;
    played(mushra_settings.expe_order(mushra_settings.expe),f+1) = 1;
else
    playNum = 1;
    nowPlaying = 0;
    played(mushra_settings.expe_order(mushra_settings.expe),1) = 1;
end

highlight_button(handles,nowPlaying);
highlight_slider(handles,nowPlaying);

finishedPlaying = play(mushra_settings.player,mushra_settings.players{mushra_settings.expe_order(mushra_settings.expe),playNum});
if finishedPlaying && ~isempty(nowPlaying)
    unhighlight_button(handles,nowPlaying);
    unhighlight_slider(handles,nowPlaying);
    nowPlaying = [];
end
%

function highlight_button(handles,button)
% eval(sprintf("set(handles.%s,'FontWeight','bold','FontSize',16);",button));
eval(sprintf("set(handles.play%d,'backgroundcolor',[0.3,1,0.3]);",button));

function unhighlight_button(handles,button)
% eval(sprintf("set(handles.%s,'FontWeight','normal','FontSize',14);",button));
eval(sprintf("set(handles.play%d,'backgroundcolor',[0.94 0.94 0.94]);",button));

function highlight_slider(handles,slider)
if slider
    eval(sprintf("set(handles.slider%d,'backgroundcolor',[0.2,0.3,1]);",slider));
end

function unhighlight_slider(handles,slider)
if slider
    eval(sprintf("set(handles.slider%d,'backgroundcolor',[0.94 0.94 0.94]);",slider));
end

function keyPressCallback_callback(hObject, eventdata, handles)
global mushra_settings nowPlaying
if contains('12345678',eventdata.Key)
    if str2double(eventdata.Key) > mushra_settings.nbfile(mushra_settings.expe_order(mushra_settings.expe))+1
        return
    end
    play_callback(handles,str2double(eventdata.Key)-1)
else
    switch eventdata.Key
        case 'q'
            stop_callback(handles)
        case 'return'
            results_callback(handles)
        case 'w'
            if nowPlaying
                move_slider(handles,nowPlaying,5);
            end
        case 'a'
            if ~isempty(nowPlaying) && nowPlaying>0
                play_callback(handles,nowPlaying-1)
            end
        case 'd'
            if ~isempty(nowPlaying) && nowPlaying < mushra_settings.nbfile(mushra_settings.expe_order(mushra_settings.expe))
                play_callback(handles,nowPlaying+1)
            end
        case 's'
            if nowPlaying
                move_slider(handles,nowPlaying,-5);
            end
    end
end

% Move slider up or down
function move_slider(handles,sliderNum,delta)
slider_obj = eval(sprintf("handles.slider%d",sliderNum));
value = get(slider_obj,'value');
value = value + delta;
if value > 100
    value = 100;
elseif value < 0
    value = 0;
end
set(slider_obj,'value',value);
slider_callback(handles,sliderNum);


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