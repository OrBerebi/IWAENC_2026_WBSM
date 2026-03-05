function plot_iwanc_results(filename)
    % plot_iwanc_results_5x2.m
    % Parses IWANC_Final_Results.txt and generates a 5x2 error-bar plot
    % Columns: Inside FoV, Outside FoV
    % Rows: stft_NMSE, stft_MagErr, stft_BSD, RMSE_BMS, RMSE_ILD
   
    FS = 22;
    if ~isfile(filename)
        error('File %s not found in the current directory.', filename);
    end
    
    %% 1. PARSE THE TEXT FILE
    fid = fopen(filename, 'r');
    data = struct();
    current_cond = '';
    algos = {};
    
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        
        % Identify Condition
        if startsWith(line, 'CONDITION:')
            current_cond = strtrim(strrep(line, 'CONDITION:', ''));
            data.(current_cond) = struct();
        end
        
        % Identify Algorithms from the Header
        if startsWith(line, 'Metric')
            parts = strsplit(line, '|');
            algos = cell(1, length(parts)-2);
            for i = 2:length(parts)-1
                algos{i-1} = strtrim(parts{i});
            end
        end
        
        % Extract Data Lines (lines containing '+/-')
        if contains(line, '+/-') && ~isempty(current_cond)
            parts = strsplit(line, '|');
            metric_name = strtrim(parts{1});
            
            for a = 1:length(algos)
                val_str = strtrim(parts{a+1});
                nums = sscanf(val_str, '%f +/- %f');
                if numel(nums) == 2
                    data.(current_cond).(metric_name).(algos{a}).mean = nums(1);
                    data.(current_cond).(metric_name).(algos{a}).std  = nums(2);
                end
            end
        end
    end
    fclose(fid);
    
    fprintf('Successfully parsed data for %d conditions.\n', length(fieldnames(data)));
    
    %% 2. PLOTTING CONFIGURATION
    fov_conditions = {'Inside', 'Outside'};
    
    % Reordered to display Inf -> 12 -> 6 -> 0
    drr_levels = {'Inf', '12dB', '6dB', '0dB'};
    drr_labels = {'\infty', '12', '6', '0'};
    
    % The 5 specific metrics (kept internal names the same for parsing)
    metrics_to_plot = {'stft_NMSE', 'stft_MagErr', 'stft_BSD', 'RMSE_BMS', 'RMSE_ILD'};
    
    % Updated the visual Y-axis labels
    y_labels = {'NMSE', 'MagErr', 'BSD', 'e_{BMS}', 'e_{ILD}'};
    
    % High-contrast, publication-safe colors
    colors = [
        0.00 0.45 0.74; % Blue (BSM)
        0.85 0.33 0.10; % Orange/Red (WBSM)
        0.47 0.67 0.19; % Green (COM)
    ];
    markers = {'s', 'v', 'o'}; % Square, Down-Triangle, Circle
    x_offsets = [-0.15, 0, 0.15]; % Stagger error bars
    
    %% 3. GENERATE THE 5x2 FIGURE
    % Set WindowStyle to docked
    fig = figure('Name', 'IWANC Final Results - 5x2 Grid', ...
                 'WindowStyle', 'docked', ...
                 'Color', 'w');
                 
    % INITIALIZE TILED LAYOUT (Replaces subplot)
    t = tiledlayout(5, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for m = 1:length(metrics_to_plot)
        metric = metrics_to_plot{m};
        
        for f = 1:length(fov_conditions)
            fov = fov_conditions{f};
            
            % Calculate subplot index (Row m, Col f)
            plot_idx = (m - 1) * 2 + f;
            
            % USE NEXTTILE INSTEAD OF SUBPLOT
            nexttile(plot_idx);
            hold on; grid on;
            set(gca, 'GridAlpha', 0.4, 'LineWidth', 1.0, 'FontSize', FS);
            
            % Loop through algorithms to plot
            for a = 1:length(algos)
                algo = algos{a};
                
                means = NaN(1, length(drr_levels));
                stds  = NaN(1, length(drr_levels));
                
                % Gather data across DRR levels
                for d = 1:length(drr_levels)
                    cond_name = sprintf('FoV_%s_DRR_%s', fov, drr_levels{d});
                    if isfield(data, cond_name) && isfield(data.(cond_name), metric)
                        % Safely extract the data
                        try
                            means(d) = data.(cond_name).(metric).(algo).mean;
                            stds(d)  = data.(cond_name).(metric).(algo).std;
                        catch
                            % Keep as NaN if algorithm is missing for this metric
                        end
                    end
                end
                
                % Plot Error Bars
                x_pos = (1:length(drr_levels)) + x_offsets(a);
                errorbar(x_pos, means, stds, ...
                    'LineStyle', 'none', ...
                    'Marker', markers{a}, ...
                    'MarkerSize', 7, ...
                    'MarkerFaceColor', 'w', ...
                    'MarkerEdgeColor', colors(a,:), ...
                    'Color', colors(a,:), ...
                    'LineWidth', 1.5, ...
                    'CapSize', 4, ...
                    'DisplayName', upper(algo));
            end
            
            % --- Subplot Formatting ---
            
            % Column headers on the top row
            if m == 1
                title(sprintf('%s FoV', fov), 'FontSize', FS, 'FontWeight', 'bold');
            end
            
            % Y-axis labels only on the left column
            if f == 1
                ylabel(y_labels{m}, 'FontSize', FS, 'FontWeight', 'bold', 'Interpreter', 'tex');
            end
            
            xlim([0.5, length(drr_levels) + 0.5]);
            
            % X-axis ticks and labels only on the bottom row
            if m == length(metrics_to_plot)
                xticks(1:length(drr_levels));
                xticklabels(drr_labels);
                set(gca, 'TickLabelInterpreter', 'tex');
                xlabel('DRR (dB)', 'FontSize', FS, 'FontWeight', 'bold');
            else
                xticks(1:length(drr_levels));
                xticklabels({});
            end
            
            % Pad Y-axis limits slightly for a cleaner look
            ylim_curr = ylim;
            y_range = ylim_curr(2) - ylim_curr(1);
            ylim([ylim_curr(1) - 0.1*y_range, ylim_curr(2) + 0.1*y_range]); 
        end
    end

    % Optional tight layout adjustment
    try
        fontsize(fig, FS, 'points'); 
    catch
    end

    %% 4. GLOBAL LEGEND AND EXPORT
    % Create a single legend that spans the entire tiled layout
    lgd = legend('Orientation', 'horizontal');
    lgd.Layout.Tile = 'south'; % This will now work perfectly!

    [folder_path, old_name, old_ext] = fileparts(filename);
    new_filename = 'results.fig';
    fig_full_path = fullfile(folder_path, new_filename);
    savefig(fig_full_path)
end