function run_beta_ablation(base_path, save_folder_path)
%RUN_BETA_ABLATION Performs ablation study on the weighting parameter beta
%   Calculates NMSE and Magnitude Error for WBSM across a range of betas.
%   Generates 2 figures (LS and MLS), each with separate NMSE and Mag subplots.

    %% 1. Setup
    % Define beta range
    % Added missing comma between 0.7 and 0.8
    beta_vals = [0.001, 0.005, 0.01, 0.05,0.075, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0];
    %beta_vals = [0.05,0.075, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0];
    %beta_vals = [0.075];
    
    n_betas = length(beta_vals);
    
    % Initialize storage vectors (Separate for LS and MLS)
    nmse_in_ls   = zeros(1, n_betas);
    nmse_out_ls  = zeros(1, n_betas);
    mag_in_ls    = zeros(1, n_betas);
    mag_out_ls   = zeros(1, n_betas);
    
    nmse_in_mls  = zeros(1, n_betas);
    nmse_out_mls = zeros(1, n_betas);
    mag_in_mls   = zeros(1, n_betas);
    mag_out_mls  = zeros(1, n_betas);

    fprintf('Starting Beta Ablation Study (%d iterations)...\n', n_betas);




    


    % array config
    arraysetup.M            = 7;
    arraysetup.arrayType    = 6;
    arraysetup.sphereType   = "rigid";
    arraysetup.fs           = 48e3; 
    [arraysetup.th_array, ~, arraysetup.ph_array, arraysetup.r_array] = ...
        BSM_toolbox.GetArrayPositions(arraysetup.arrayType, arraysetup.M, 0);

    fprintf('Pre-loading data structures...\n');
    dummy_params = init_params_local(1, base_path, save_folder_path);
    dummy_params.fs = arraysetup.fs;
    data_static      = load_and_prepare_data_local(dummy_params); 

        % Frequency cutoff logic for MLS
    f_vec = linspace(0, round(data_static.fs/2), round(data_static.nFFT/2)+1);
    f_cutoff = 1.5e3;

    f_cutoff_high = 8e3;

    [~,idx_cutoff] = min(abs(f_vec - f_cutoff));
    [~,idx_cutoff_high] = min(abs(f_vec - f_cutoff_high));


    %% 3. Main Ablation Loop
    for i = 1:n_betas
        current_beta = beta_vals(i);
        fprintf('Processing Beta = %.4f (%d/%d)...\n', current_beta, i, n_betas);
        
        % Run the analysis (using v4 as requested)
        % Ensure anechoic_analysis_FOV_BSM_v4 accepts (beta, base, save, data, plot_flag)
        plot_anechoic_analysis = false;
        %res = anechoic_analysis_FOV_BSM_v4(current_beta, base_path, save_folder_path, data_static, false);
        res  = anechoic_analysis_FOV_BSM_v5(current_beta,arraysetup, base_path, save_folder_path, data_static, plot_anechoic_analysis);
        
        % --- Extract Broadband Errors ---
        
        % 1. MLS Analysis (Split band)
        % Low Freq for NMSE
        %nmse_in_mls(i)  = mean(mean(res.errors.nmse.in.mls.WBSM.mu(1:idx_cutoff,:)));
        %nmse_out_mls(i) = mean(mean(res.errors.nmse.out.mls.WBSM.mu(1:idx_cutoff,:)));
        
        % High Freq for Mag Error
        %mag_in_mls(i)  = mean(mean(res.errors.mag.in.mls.WBSM.mu(idx_cutoff:end,:)));
        %mag_out_mls(i)  = mean(mean(res.errors.mag.out.mls.WBSM.mu(idx_cutoff:end,:)));

        nmse_in_mls(i)  = mean(mean(10*log10(res.errors.nmse.in.mls.WBSM.mu(1:idx_cutoff,:)),2),1);
        nmse_out_mls(i) = mean(mean(10*log10(res.errors.nmse.out.mls.WBSM.mu(1:idx_cutoff,:)),2),1);
        mag_in_mls(i)   = mean(mean(10*log10(res.errors.mag.in.mls.WBSM.mu(idx_cutoff:idx_cutoff_high,:)),2),1);
        mag_out_mls(i)  = mean(mean(10*log10(res.errors.mag.out.mls.WBSM.mu(idx_cutoff:idx_cutoff_high,:)),2),1);
        
        % 2. LS Analysis (Full band)
        % nmse_in_ls(i)   = mean(mean(res.errors.nmse.in.ls.WBSM.mu));
        % nmse_out_ls(i)  = mean(mean(res.errors.nmse.out.ls.WBSM.mu));
        % mag_in_ls(i)    = mean(mean(res.errors.mag.in.ls.WBSM.mu));
        % mag_out_ls(i)   = mean(mean(res.errors.mag.out.ls.WBSM.mu));

        nmse_in_ls(i)  = mean(mean(10*log10(res.errors.nmse.in.ls.WBSM.mu),2),1);
        nmse_out_ls(i) = mean(mean(10*log10(res.errors.nmse.out.ls.WBSM.mu),2),1);
        mag_in_ls(i)   = mean(mean(10*log10(res.errors.mag.in.ls.WBSM.mu),2),1);
        mag_out_ls(i)  = mean(mean(10*log10(res.errors.mag.out.ls.WBSM.mu),2),1);

    end

    %% 4. Plotting
    
    ablation_save_path = fullfile(base_path, save_folder_path, "Ablation_Study");
    if ~exist(ablation_save_path, 'dir')
        mkdir(ablation_save_path);
    end

    % --- Figure 1: LS Solver (Top: NMSE, Bottom: Mag) ---
    fig_ls = figure('Name', 'Ablation: LS Solver', 'Color', 'w');
    t_ls = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Top Tile: NMSE
    nexttile;
    semilogx(beta_vals, nmse_in_ls, '-o', 'LineWidth', 2, 'DisplayName', 'Within FoV');
    hold on;
    semilogx(beta_vals, nmse_out_ls, '-s', 'LineWidth', 2, 'DisplayName', 'Outside FoV');
    grid on;
    ylabel('NMSE [dB]');
    title('LS: Broadband NMSE');
    legend('Location', 'best');
    xlim([min(beta_vals), max(beta_vals)]);
    
    % Bottom Tile: Magnitude Error
    nexttile;
    semilogx(beta_vals, mag_in_ls, '-o', 'LineWidth', 2, 'DisplayName', 'Within FoV');
    hold on;
    semilogx(beta_vals, mag_out_ls, '-s', 'LineWidth', 2, 'DisplayName', 'Outside FoV');
    grid on;
    xlabel('Beta (\beta)');
    ylabel('Mag Error [dB]');
    title('LS: Broadband Magnitude Error');
    % legend('Location', 'best'); % Optional to repeat legend
    xlim([min(beta_vals), max(beta_vals)]);
    xlabel('$$\beta$$');
    
    
    sgtitle('Ablation Study: LS Solver');
    
    % Save LS Figure
    saveas(fig_ls, fullfile(ablation_save_path, 'Ablation_LS_Combined.fig'));
    saveas(fig_ls, fullfile(ablation_save_path, 'Ablation_LS_Combined.png'));

    
    % --- Figure 2: MLS Solver (Top: NMSE, Bottom: Mag) ---
    fig_mls = figure('Name', 'Ablation: MLS Solver', 'Color', 'w', 'Position', [750 100 600 800]);
    t_mls = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Top Tile: NMSE (Low Freq)
    nexttile;
    semilogx(beta_vals, nmse_in_mls, '-o', 'LineWidth', 2, 'DisplayName', 'Within FoV');
    hold on;
    semilogx(beta_vals, nmse_out_mls, '-s', 'LineWidth', 2, 'DisplayName', 'Outside FoV');
    grid on;
    ylabel('NMSE [dB]');
    title(sprintf('MLS: NMSE (f bellow %.1f kHz)', f_cutoff/1e3));
    legend('Location', 'best');
    xlim([min(beta_vals), max(beta_vals)]);
    
    % Bottom Tile: Magnitude Error (High Freq)
    nexttile;
    semilogx(beta_vals, mag_in_mls, '-o', 'LineWidth', 2, 'DisplayName', 'Within FoV');
    hold on;
    semilogx(beta_vals, mag_out_mls, '-s', 'LineWidth', 2, 'DisplayName', 'Outside FoV');
    grid on;
    xlabel('Beta (\beta)');
    ylabel('Mag Error [dB]');
    title(sprintf('MLS: Magnitude Error (f above %.1f kHz)', f_cutoff/1e3));
    xlim([min(beta_vals), max(beta_vals)]);
    xlabel('$$\beta$$');
    
    sgtitle('Ablation Study: MLS Solver');
    
    % Save MLS Figure
    saveas(fig_mls, fullfile(ablation_save_path, 'Ablation_MLS_Combined.fig'));
    saveas(fig_mls, fullfile(ablation_save_path, 'Ablation_MLS_Combined.png'));
    
    fprintf('Ablation study complete. Results saved to %s\n', ablation_save_path);
end

%% ================= COPIED HELPERS FOR DATA LOADING =================
function params = init_params_local(beta, base_path, save_folder_path)
    params.beta   = beta;
    params.base   = base_path;
    params.save   = base_path + save_folder_path + "/beta_" + num2str(beta) + "/";
    params.HRTFpath    = base_path + "/HRTF/HRIR_L2702_koln.sofa";
    params.ATFpath     = base_path + "/ATF/aria_interpolated.mat";
end

function data = load_and_prepare_data_local(params)
    startup_script(params.base)
    hobj       = sofaToEaro(convertStringsToChars(params.HRTFpath));
    omega_hrtf = [hobj.sourceGrid.elevation,hobj.sourceGrid.azimuth,hobj.sourceGrid.r];

    if hobj.fs ~= params.fs
        hobj = hobj.resampleData(params.fs);
        %hobj.fs = params.fs;
    end
    ATF = load(params.ATFpath).aria;
    ATF.IR = ATF.IR_interpolated;
    %ATF.IR = permute(ATF.IR,[1,3,2]);

     if ATF.fs ~= params.fs
        gComDiv = gcd(double(ATF.fs), params.fs);
        p = double(params.fs / gComDiv);
        q = double(ATF.fs / gComDiv);
        newData = zeros(size(ATF.IR,1), ceil(size(ATF.IR,2) * (p/q)), size(ATF.IR,3));
         % Resample the SH RIRs (Optional, but good for completeness)
        for ii=1:size(ATF.IR,1)
           newData(ii,:,:) = resample(squeeze(ATF.IR(ii,:,:)),p,q);
        end
        ATF.IR = newData;

    end
    ATF.IR = permute(ATF.IR,[1,3,2]);
    az_res = 1;
    phi_az = deg2rad(linspace(0,180,180/az_res).');
    th_az  = deg2rad(90*ones(size(phi_az)));
    r_az   = ones(size(phi_az));
    omega_az = [th_az, phi_az, r_az];

    tmp1 = 2^(nextpow2(size(ATF.IR,3)));
    tmp2 = 2^(nextpow2(size(hobj.data,2)));
    nFFT = max([tmp1,tmp2,2^10]);
    %nFFT = 2^10;

    ATF.IR = cat(3, ATF.IR, zeros(size(ATF.IR,1),size(ATF.IR,2),nFFT-size(ATF.IR,3)));
    hobj.data = cat(2,hobj.data,zeros(size(hobj.data,1),nFFT-size(hobj.data,2),size(hobj.data,3)));
    hobj.taps = nFFT;

    data.hobj = hobj;
    data.ATF = ATF;
    data.omega_hrtf = omega_hrtf;
    data.omega_az   = omega_az;
    data.nFFT = nFFT;
    data.fs   = hobj.fs;
end