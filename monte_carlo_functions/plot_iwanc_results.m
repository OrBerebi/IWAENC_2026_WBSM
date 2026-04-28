function plot_iwanc_results(filename)
    % plot_iwanc_results_6x2.m
    % Parses IWANC_Final_Results.txt and generates a 6x2 error-bar plot
   
    FS = 22;
    if ~isfile(filename)
        error('File %s not found in the current directory.', filename);
    end
    
    %% 1. PARSE THE TEXT FILE (Unchanged)
    fid = fopen(filename, 'r');
    data = struct();
    current_cond = '';
    algos = {};
    
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if startsWith(line, 'CONDITION:')
            current_cond = strtrim(strrep(line, 'CONDITION:', ''));
            data.(current_cond) = struct();
        end
        if startsWith(line, 'Metric')
            parts = strsplit(line, '|');
            algos = cell(1, length(parts)-2);
            for i = 2:length(parts)-1
                algos{i-1} = strtrim(parts{i});
            end
        end
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
    
    %% 2. PLOTTING CONFIGURATION
    fov_conditions = {'Inside', 'Outside'};
    drr_levels = {'Inf', '12dB', '6dB', '0dB'};
    drr_labels = {'\infty', '12', '6', '0'};
    
    % Metrics and their corresponding Y-labels
    metrics_to_plot = {'stft_NMSE', 'stft_MagErr', 'stft_BSD', 'RMSE_BMS', 'RMSE_ILD','RMSE_IC'};
    y_labels = {'NMSE', 'MagErr', 'BSD', 'e_{BMS}', 'e_{ILD}','e_{IC}'};
    
    % --- NEW: DEFINED Y-LIMITS ---
    ylim_values = {
        [-40, 0], ... % NMSE
        [-25, 0], ... % MagErr
        [0, 8],   ... % BSD
        [0, 8],   ... % e_BMS
        [0, 8],   ... % e_ILD
        [0, 1]    ... % e_IC
    };
    
    colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];
    markers = {'s', 'v', 'o'};
    x_offsets = [-0.15, 0, 0.15];
    
    %% 3. GENERATE THE 6x2 FIGURE
    fig = figure('Name', 'IWANC Final Results', 'WindowStyle', 'docked', 'Color', 'w');
    t = tiledlayout(6, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for m = 1:length(metrics_to_plot)
        metric = metrics_to_plot{m};
        
        for f = 1:length(fov_conditions)
            fov = fov_conditions{f};
            plot_idx = (m - 1) * 2 + f;
            
            nexttile(plot_idx);
            hold on; grid on;
            set(gca, 'GridAlpha', 0.4, 'LineWidth', 1.0, 'FontSize', FS);
            
            for a = 1:length(algos)
                algo = algos{a};
                means = NaN(1, length(drr_levels));
                stds  = NaN(1, length(drr_levels));
                
                for d = 1:length(drr_levels)
                    cond_name = sprintf('FoV_%s_DRR_%s', fov, drr_levels{d});
                    if isfield(data, cond_name) && isfield(data.(cond_name), metric)
                        try
                            means(d) = data.(cond_name).(metric).(algo).mean;
                            stds(d)  = data.(cond_name).(metric).(algo).std;
                        catch
                        end
                    end
                end
                
                x_pos = (1:length(drr_levels)) + x_offsets(a);
                errorbar(x_pos, means, stds, 'LineStyle', 'none', ...
                    'Marker', markers{a}, 'MarkerSize', 7, 'MarkerFaceColor', 'w', ...
                    'MarkerEdgeColor', colors(a,:), 'Color', colors(a,:), ...
                    'LineWidth', 1.5, 'CapSize', 4, 'DisplayName', upper(algo));
            end
            
            % --- Updated Formatting with Constant Y-Limits ---
            if m == 1
                title(sprintf('%s FoV', fov), 'FontSize', FS, 'FontWeight', 'bold');
            end
            
            if f == 1
                ylabel(y_labels{m}, 'FontSize', FS, 'FontWeight', 'bold', 'Interpreter', 'tex');
            end
            
            xlim([0.5, length(drr_levels) + 0.5]);
            
            % Apply the specific constant ylim for this row
            ylim(ylim_values{m}); 
            
            if m == length(metrics_to_plot)
                xticks(1:length(drr_levels));
                xticklabels(drr_labels);
                xlabel('DRR (dB)', 'FontSize', FS, 'FontWeight', 'bold');
            else
                xticks(1:length(drr_levels));
                xticklabels({});
            end
        end
    end

    %% 4. GLOBAL LEGEND AND EXPORT
    lgd = legend('Orientation', 'horizontal');
    lgd.Layout.Tile = 'south'; 
    [folder_path, ~, ~] = fileparts(filename);
    savefig(fullfile(folder_path, 'results.fig'));
end