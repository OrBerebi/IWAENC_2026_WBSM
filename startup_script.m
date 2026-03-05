function startup_script(base_path)

    
    % Set graphics defaults:
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot,'defaulttextinterpreter','latex');
    set(0,'defaultTextInterpreter','latex'); % Keep default as LaTeX interpreter
    % Set line properties:
    set(0, 'DefaultLineLineWidth', 2);
    set(0, 'DefaultLineMarkerSize', 10);
    % Enable grids:
    set(0, 'defaultAxesXGrid', 'on');
    set(0, 'defaultAxesYGrid', 'on');
    set(0, 'defaultAxesZGrid', 'on');
    % Enable minor grids:
    set(groot, 'defaultAxesXMinorGrid', 'on', 'defaultAxesXMinorGridMode', 'manual');
    set(groot, 'defaultAxesYMinorGrid', 'on', 'defaultAxesYMinorGridMode', 'manual');
    set(groot, 'defaultAxesZMinorGrid', 'on', 'defaultAxesZMinorGridMode', 'manual');
    % Set font properties:
    set(0, 'DefaultTextFontName', 'Times New Roman');
    set(0, 'DefaultAxesFontName', 'Times New Roman');
    set(0, 'DefaultTextFontSize', 18); % Sets font size to 
    set(0, 'defaultAxesFontSize', 18); % Sets axes font size to
    set(0, 'DefaultAxesTitleFontSizeMultiplier', 1); % Keep title font size multiplier
    % Set figure window style:
    set(0, 'DefaultFigureWindowStyle', 'docked');
    
    addpath(genpath(base_path));
end