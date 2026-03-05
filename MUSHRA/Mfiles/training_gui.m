function fig=training_gui(nbfile,expe,nbexpe,initPos, metricString)

% EVALUATION_GUI Creates the GUI for the evaluation phase
%
% fig=evaluation_gui(nbfile,expe,nbexpe) creates the GUI given the number
% of tested files nbfile (including the reference), the experiment number
% expe and the total number of experiments nbexpe
%
% returns a handle to the figure if succeeds (nbfile small enough),
% returns false otherwise
UseinitPos = true;
if ~exist('initPos','var') || isempty(initPos)
    UseinitPos = false;
end
if ~exist('metricString','var') || isempty(metricString)
    metricString = '';
end

%%%getting the screen size and computing the figure width
%screen size
set(0,'Units','characters');
siz=get(0,'Screensize');
maxwidth=siz(3);
FS = 0.6;
%width of each item
iwidth=min(14,(maxwidth-26.5)/nbfile);
if iwidth < 9.5
    fig=[];
    return;
end
ihight = 1.5;
%figure width
width=50+iwidth*nbfile; %37

%%%opening the figure
if UseinitPos
    PositionFig = [initPos width 32];
else
    PositionFig = [0 0 width 32];
end
fig=figure('Name','MUSHRAM - Training phase','NumberTitle','off','MenuBar','none','Resize','off','Units','characters','Position',PositionFig);
movegui(fig,'center')
%%%displaying fixed items
%title
uicontrol(fig,'Style','Text','Units','characters','Position',[0 29 width ihight],'FontUnits','normalized','FontSize',FS,'String',['Training ' num2str(expe) '/' int2str(nbexpe)],'Tag','experiment');
uicontrol(fig,'Style','Text','Units','characters','Position',[0 27.5 width ihight],'FontUnits','normalized','FontSize',FS,'String',['Rate the similarity between the test signals '],'Tag','experiment2');
uicontrol(fig,'Style','Text','Units','characters','Position',[0 26.2 width ihight],'FontUnits','normalized','FontSize',FS,'String',['and the reference with respect to ' metricString '.'],'Tag','experiment3');

%verbal performance assessments
uicontrol(fig,'Style','Text','Units','characters','Position',[3.5 21.57 19 ihight],'FontUnits','normalized','FontSize',FS,'String','Excellent','Tag','scale90');
uicontrol(fig,'Style','Text','Units','characters','Position',[3.5 19.24 19 ihight],'FontUnits','normalized','FontSize',FS,'String','Good','Tag','scale70');
uicontrol(fig,'Style','Text','Units','characters','Position',[3.5 16.92 19 ihight],'FontUnits','normalized','FontSize',FS,'String','Fair','Tag','scale50');
uicontrol(fig,'Style','Text','Units','characters','Position',[3.5 14.59 19 ihight],'FontUnits','normalized','FontSize',FS,'String','Poor','Tag','scale30');
uicontrol(fig,'Style','Text','Units','characters','Position',[3.5 12.27 19 ihight],'FontUnits','normalized','FontSize',FS,'String','Bad','Tag','scale10');
%save and proceed button
results=uicontrol(fig,'Style','Pushbutton','Units','characters','FontUnits','normalized','FontSize',FS-0.15,'Tag','results','Callback','training_callbacks(''results'',guidata(gcbo))');
set(results,'Position',[width/2-8 1.5 25 2],'String','Save & proceed (enter)');

%play reference button
uicontrol(fig,'Style','Pushbutton','Units','characters','Position',[3.5 5.7 19 2],'FontUnits','normalized','FontSize',FS-0.15,'String','Play reference (1)','Callback','training_callbacks(''play'',guidata(gcbo),0)','Tag','play0');
%horizontal dashed lines
ax=axes('Units','characters','Position',[-.5 -.5 width+.5 32.5],'Tag','axes');
plot(ax,[3.5 width-3.5],[12.15 12.15],'--k',[3.5 width-3.5],[14.51 14.51],'--k',[3.5 width-3.5],[16.87 16.87],'--k',[3.5 width-3.5],[19.23 19.23],'--k',[3.5 width-3.5],[21.59 21.59],'--k',[3.5 width-3.5],[23.95 23.95],'--k');
set(ax,'Color',[0.701961 0.701961 0.701961],'XTick',[],'YTick',[],'XLim',[-.5 width+.5],'YLim',[-.5 32.5]);

% Stop button
uicontrol(fig,'Style','Pushbutton','Units','characters','Position',[3.5 8 9 2],'FontUnits','normalized','FontSize',FS-0.15,'String','Stop (q)','Tag','Stop','Callback','training_callbacks(''stop'',guidata(gcbo))');

%%%displaying items depending on the number of test signals
for f=1:nbfile
    %sound number
    uicontrol(fig,'Style','Text','Units','characters','Position',[14+f*iwidth 24.5 9 ihight],'FontUnits','normalized','FontSize',FS,'String',['Sound ' char(64+f)],'Tag',['title' int2str(f)]);
    %slider
    uicontrol(fig,'Style','Slider','Units','characters','Position',[17.1+f*iwidth 11.24 2.8 13],'Min',0,'Max',100,'Value',0,'Callback',['training_callbacks(''slider'',guidata(gcbo),' int2str(f) ')'],'Tag',['slider' int2str(f)]);
    %performance figures
    uicontrol(fig,'Style','Edit','Enable','inactive','Units','characters','Position',[14+f*iwidth 8.2 9 1.6],'FontUnits','normalized','FontSize',FS,'String',00,'Tag',['rating' int2str(f)]);
    %play buttons
    uicontrol(fig,'Style','Pushbutton','Units','characters','Position',[14+f*iwidth 5.7 9 2],'FontUnits','normalized','FontSize',FS-0.15,'String',['Play (' int2str(f+1) ')'],'Callback',['training_callbacks(''play'',guidata(gcbo),' int2str(f) ')'],'Tag',['play' int2str(f)]);
end

% textbox for keyboard shortcuts
keyTitle = uicontrol(fig,'Style','Text','Units','characters','Position',[0 0.3 26 5],'FontUnits','normalized','FontSize',FS-0.45,'String',sprintf('Keyboard Shortcuts'),'Tag','keyboardShortcut');
key = uicontrol(fig,'Style','Text','Units','characters','Position',[1 0 26 4.1],'FontUnits','normalized','FontSize',FS-0.45,'String',sprintf([' Play  (#1-#' num2str(nbfile+1) ')    Stop             (q)\n Move left (a)    Move right     (d)\n Slider up  (w)    Slider down   (s)']),'HorizontalAlignment','left','Tag','keyboardShortcut');
set(key,'BackgroundColor',[0.701961 0.701961 0.701961])
set(keyTitle,'BackgroundColor',[0.701961 0.701961 0.701961])
ax=axes('Units','characters','Position',[0 0 27 5.5],'Tag','axes2');
plot(ax,[0 26],[5 5],'-k',[26 26],[0 5],'-k');
set(ax,'Color',[0.701961 0.701961 0.701961],'XTick',[],'YTick',[]);

set(fig,'WindowKeyPressFcn',@(hObject,eventdata)training_callbacks('keyPressCallback',hObject,eventdata,guidata(hObject)));
set(fig,'CloseRequestFcn',@(hObject,eventdata)training_callbacks('CloseFigure',hObject,eventdata,guidata(hObject)));

%%%moving onscreen
movegui(fig,'onscreen');
