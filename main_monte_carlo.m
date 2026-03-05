% main_run_iwanc_simulations.m
clear; clc; 
close all;

base_path = string(pwd);
startup_script(base_path); 
results_dir_path = base_path + "/results/04_March_2026_Aria_M5/";
mkdir(results_dir_path);
results_txt_file = fullfile(results_dir_path, 'IWANC_Final_Results.txt');

% Clear previous results file if starting fresh
if isfile(results_txt_file)
    delete(results_txt_file);
end

% Define Conditions
fov_conditions = ["Inside", "Outside"];

% Use 999 to represent 'Infinity' (Direct-only)
drr_conditions = [999, 12, 6, 0]; 
%drr_conditions = [0]; 
room_types = ["small","medium","large"];
target_counts = [1,2,3];

fprintf('Starting IWANC Automated Pipeline...\n');

for f = 1:length(fov_conditions)
    for d = 1:length(drr_conditions)
        
        fov_cond = fov_conditions(f);
        drr_val = drr_conditions(d);
        
        if drr_val == 999
            cond_name = sprintf('FoV_%s_DRR_Inf', fov_cond);
        else
            cond_name = sprintf('FoV_%s_DRR_%ddB', fov_cond, drr_val);
        end
        
        fprintf('\n******************************************************\n');
        fprintf('RUNNING CONDITION: %s\n', cond_name);
        fprintf('******************************************************\n');
        
        % Define isolated directories for this specific condition
        data_dir       = fullfile(results_dir_path, sprintf('simulation_data_%s', cond_name));
        results_dir    = fullfile(results_dir_path, sprintf('evaluation_results_%s', cond_name));
        scenarios_file = fullfile(results_dir_path, sprintf('simulation_scenarios_%s.mat', cond_name));
        
        % --- STAGE 1: Generate ---
        monte_carlo_01_generate_data(results_dir_path, drr_val, fov_cond, room_types ,target_counts, data_dir, scenarios_file);
        
        % --- STAGE 2: Evaluate ---
        monte_carlo_02_evaluate_data(base_path, data_dir, scenarios_file, results_dir);
        
        % --- STAGE 3: Analyze ---
        % Set plot_flag to false to prevent hundreds of windows from popping up
        plot_flag = false; 
        monte_carlo_03_analysis_data(base_path, results_dir, results_txt_file, cond_name, plot_flag);
        
    end
end

fprintf('\n=== ALL SIMULATIONS COMPLETE ===\n');
fprintf('Results saved to: %s\n', results_txt_file);
plot_iwanc_results(results_txt_file)
