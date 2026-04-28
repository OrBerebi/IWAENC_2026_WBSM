% =========================================================================
% True-to-Paper Figure Formatter
% Fixes Docked Window errors, matches LaTeX rendered text size, and silences warnings
% =========================================================================
clear; clc; close all;

% --- Paper Specifications (IEEE/AVAR Template) ---
FONT_NAME = 'Times New Roman';
FONT_SIZE = 20; % Standard main text size in the paper

% Standard column widths in inches
SINGLE_COLUMN_WIDTH = 3.5;  
DOUBLE_COLUMN_WIDTH = 7.16; 

% Find all .fig files in the current directory
figFiles = dir('*.fig');

if isempty(figFiles)
    disp('No .fig files found in the current directory!');
    return;
end

disp(['Found ', num2str(length(figFiles)), ' figures to process...']);

for i = 1:length(figFiles)
    filename = figFiles(i).name;
    disp(['Processing: ', filename]);
    
    % --- CRITICAL FIX: Silence MATLAB's internal interpreter warnings ---
    warning('off', 'all'); % Turn off all warnings temporarily
    fig = openfig(filename);
    warning('on', 'all');  % Turn them back on immediately after loading
    
    % Undock the figure so it can be resized
    set(fig, 'WindowStyle', 'normal');
    
    % Force figure background to white
    set(fig, 'Color', 'w');
    
    % --- Smart Figure Resizing for Pixel Perfection ---
    % Read current dimensions to figure out if it's a wide plot
    set(fig, 'Units', 'inches');
    pos = fig.Position;
    aspectRatio = pos(3) / pos(4); % Width / Height
    
    if aspectRatio > 1.6
        % It's a wide plot (like your 6-panel grid), scale to double-column
        newWidth = DOUBLE_COLUMN_WIDTH;
    else
        % It's a standard plot (like the bar chart), scale to single-column
        newWidth = SINGLE_COLUMN_WIDTH;
    end
    
    % Calculate new height while preserving the original aspect ratio
    newHeight = newWidth / aspectRatio;
    
    % Apply the new physical size (This dictates the final pixel size at 300 DPI)
    fig.Position = [pos(1), pos(2), newWidth, newHeight];
    
    % --- Update Fonts Globally ---
    % Apply 10pt Times New Roman to absolutely every text object
    textObjects = findall(fig, '-property', 'FontName');
    set(textObjects, 'FontName', FONT_NAME);
    
    sizeObjects = findall(fig, '-property', 'FontSize');
    set(sizeObjects, 'FontSize', FONT_SIZE);
    
    % --- Save & Export ---
    [~, name, ~] = fileparts(filename);
    newFigName = [name, '_PaperReady.fig'];
    exportName = [name, '_PaperReady.png'];
    
    % Save updated .fig safely
    savefig(fig, newFigName);
    
    % Export to high-res 300 DPI image
    exportgraphics(fig, exportName, 'Resolution', 300);
    
    % Close figure to free memory
    close(fig);
    disp(['   -> Saved as ', exportName]);
end

disp('All figures successfully updated and sized for LaTeX!');