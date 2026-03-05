function fig=training2_gui(nbfile,nbexpe,run_all,initPos, metricString)

% TRAINING_GUI Creates the GUI for the training phase
%
% fig=training_gui(nbfile,nbexpe,run_all) creates the GUI given the number
% of tested files nbfile (including the reference), the number of
% experiments nbexpe and the flag run_all describing whether the evaluation
% phase is to be run afterwards
%
% returns a handle to the figure if succeeds (nbfile and nbexpe small
% enough), returns false otherwise
UseinitPos = true;
if ~exist('initPos','var') || isempty(initPos)
    UseinitPos = false;
end
if ~exist('metricString','var') || isempty(metricString)
    metricString = '';
end

%%%getting the screen size and computing the figure width/height
%screen size
set(0,'Units','characters');
siz=get(0,'Screensize');
maxwidth=siz(3);
maxheight=siz(4)*.92;
%width of each item
iwidth=min(11,(maxwidth-48)/(max(nbfile)-1));
if iwidth < 9.5
    fig=false;
    return;
end
FS = 0.6;
%figure width
width=48+iwidth*(max(nbfile)-1);
%height of each item
iheight=min(2.8,(maxheight-12.3)/nbexpe);
if iheight < 2
    fig=false;
    return;
end
%figure height
height=12.3+iheight*nbexpe;

%%%opening the figure
fig=figure('Name','MUSHRAM - Signals familiarization phase','NumberTitle','off','MenuBar','none','Resize','off','Color',[0.701961 0.701961 0.701961],'Units','characters','Position',[0 0 width height]);
movegui(fig,'center')
%%%displaying fixed items
%titles
uicontrol(fig,'Style','Text','Units','characters','Position',[0 9.3+nbexpe*iheight width 1.5],'FontUnits','normalized','FontSize',FS+0.05,'String','Signals Familiarization Phase','Tag','training');
uicontrol(fig,'Style','Text','Units','characters','Position',[26.5 6.3+nbexpe*iheight 11 1.5],'FontUnits','normalized','FontSize',FS,'String','Reference','Tag','reference');
uicontrol(fig,'Style','Text','Units','characters','Position',[46.5 6.3+nbexpe*iheight width-49 1.5],'FontUnits','normalized','FontSize',FS,'String','Test','Tag','test');
%exit or proceed button
proceed=uicontrol(fig,'Style','Pushbutton','Units','characters','FontUnits','normalized','FontSize',FS-0.15,'Tag','proceed','Callback','training2_callbacks(''proceed'',guidata(gcbo))');
if run_all
    set(proceed,'Position',[width/2-3 1.5 29 2],'String','Proceed to evaluation (enter)');
else
    set(proceed,'Position',[width/2-6 1.5 12 2],'String','Exit (enter)');
end
% Stop button
uicontrol(fig,'Style','Pushbutton','Units','characters','Position',[width/2-13 1.5 9 2],'FontUnits','normalized','FontSize',FS-0.15,'String','Stop (q)','Tag','Stop','Callback','training2_callbacks(''stop'',guidata(gcbo))');

%%%displaying items depending on the number of experiments and test signals
for e=1:nbexpe
    %experiment number
    uicontrol(fig,'Style','Text','Units','characters','Position',[1.5 6.1+(nbexpe-e)*iheight 15.5 1.5],'FontUnits','normalized','FontSize',FS,'String',['Experiment ' int2str(e)],'Tag',['title' int2str(e)]);
    %reference play buttons
    uicontrol(fig,'Style','Pushbutton','Units','characters','Position',[22.5 5.7+(nbexpe-e)*iheight 19 2],'FontUnits','normalized','FontSize',FS-0.15,'String','Play reference','Callback',['training2_callbacks(''play'',guidata(gcbo),' int2str(e) ',0)'],'Tag',['play' int2str(e) '_0']);
    for f=1:nbfile(e)-1
        %play buttons
        uicontrol(fig,'Style','Pushbutton','Units','characters','Position',[46.5+(f-1)*iwidth 5.7+(nbexpe-e)*iheight 9 2],'FontUnits','normalized','FontSize',FS-0.15,'String','Play','Callback',['training2_callbacks(''play'',guidata(gcbo),' int2str(e) ',' int2str(f) ')'],'Tag',['play' int2str(e) '_' int2str(f)]);
    end
end
set(fig,'WindowKeyPressFcn',@(hObject,eventdata)training2_callbacks('keyPressCallback',hObject,eventdata,guidata(hObject)));
set(fig,'CloseRequestFcn',@(hObject,eventdata)training2_callbacks('CloseFigure',hObject,eventdata,guidata(hObject)));

% textbox for keyboard shortcuts
keyTitle = uicontrol(fig,'Style','Text','Units','characters','Position',[0 0.3 26 5],'FontUnits','normalized','FontSize',FS-0.45,'String',sprintf('Keyboard Shortcuts'),'Tag','keyboardShortcut');
key = uicontrol(fig,'Style','Text','Units','characters','Position',[1 0 26 4.1],'FontUnits','normalized','FontSize',FS-0.45,'String',sprintf([' Stop        (q)\n Move left (a)   Move right  (d)\n Move up  (w)  Move down (s)']),'HorizontalAlignment','left','Tag','keyboardShortcut');
set(key,'BackgroundColor',[0.701961 0.701961 0.701961])
set(keyTitle,'BackgroundColor',[0.701961 0.701961 0.701961])
ax=axes('Units','characters','Position',[0 0 27 5.5],'Tag','axes2');
plot(ax,[0 26],[5 5],'-k',[26 26],[0 5],'-k');
set(ax,'Color',[0.701961 0.701961 0.701961],'XTick',[],'YTick',[]);


%%%moving onscreen
movegui(fig,'onscreen');