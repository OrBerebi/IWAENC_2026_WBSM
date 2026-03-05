function res_out = anechoic_analysis_FOV_BSM_v5(beta,arraysetup, base_path, save_folder_path, data_in, plot_flag)
%ANECHOIC_ANALYSIS_FOV_BSM_V2 Main pipeline for BSM analysis
%   Modular version of the original long script

    if nargin < 4, data_in = []; end
    if nargin < 5, plot_flag = false; end

    %% --- Parameters & Paths ---
    params = init_params(beta,arraysetup, base_path, save_folder_path);
    params.plot_figures = plot_flag; % OVERRIDE the internal true setting

    %% --- Data Loading ---
    if isempty(data_in)
        data = load_and_prepare_data(params);
    else
        data = data_in;
    end

    %% --- Setup BSM Objects ---
    [objs,data] = setup_BSM_objects(params, data);

    %% --- Compute Filters ---
    [coeffs,data] = compute_BSM_filters(objs, data, params);


    hobj_f_lebedev      = data.hobj;
    hobj_f_lebedev.data = circshift(hobj_f_lebedev.data,hobj_f_lebedev.taps/2 ,2);
    hobj_f_lebedev      = hobj_f_lebedev.toSH(44,'SRC');
    hobj_f_lebedev      = hobj_f_lebedev.toSpace('SRC',data.omega_bsm(:,1),data.omega_bsm(:,2));
    hobj_f_lebedev.data = real(hobj_f_lebedev.data);
    hobj_f_lebedev      = hobj_f_lebedev.toFreq(data.nFFT);
    hobj_f_lebedev.data = hobj_f_lebedev.data(:, 1:ceil(data.nFFT/2)+1, :);
    HRTF = hobj_f_lebedev.data;

    lebedev_order = 35;
    [gridData_Inf, ~, ~] = sofia_lebedev(Order_to_degree_lebedev(lebedev_order),0);

    lebedev_weights = gridData_Inf(:,3);
    % coeffs.BSM.f_ls  = apply_covariance_constraint_lebedev(coeffs.BSM.f_ls, data.V_k_bsm, HRTF,lebedev_weights);
    % coeffs.BSM.f_mls = apply_covariance_constraint_lebedev(coeffs.BSM.f_mls, data.V_k_bsm, HRTF,lebedev_weights);
    % coeffs.TBSM.f_ls  = apply_covariance_constraint_lebedev(coeffs.TBSM.f_ls, data.V_k_bsm, HRTF,lebedev_weights);
    % coeffs.TBSM.f_mls = apply_covariance_constraint_lebedev(coeffs.TBSM.f_mls, data.V_k_bsm, HRTF,lebedev_weights);

    % Vf = data.V_k_bsm;
    % Hf = HRTF;
    % f_sig = linspace(0,data.fs/2,data.nFFT/2+1);
    % f_sig(1) = f_sig(2)/4;
    % algos_to_plot = ["BSM","TBSM"];
    % plot_spatial_diffuse_IC(coeffs, Vf, Hf, f_sig, lebedev_weights, algos_to_plot)


    %% --- Estimate HRTFs ---
    hrtf_hat = estimate_hrtf(objs, coeffs, data, params);

    %% --- Evaluate Errors ---
    [errors,data] = evaluate_errors(hrtf_hat, data, objs, params);
    
    %% --- Plot Results ---
    if params.plot_figures
        mkdir(params.save);
        plot_results(errors, hrtf_hat, data, objs, params);
    end
    %% --- Package Output ---
    res_out = package_results(objs, data, coeffs, errors,hrtf_hat);

end

%% ===================== HELPER FUNCTIONS =====================


function plot_spatial_diffuse_IC(coeffs, Vf, Hf, f_sig, lebedev_weights, algos_to_plot)
    % PLOT_SPATIAL_DIFFUSE_IC 
    % Calculates the Interaural Coherence (IC) of the spatial field by 
    % integrating across the 1730 Lebedev directions.
    %
    % Inputs:
    %   coeffs          - Struct containing BSM weights (e.g., coeffs.BSM) [Q x 2 x nBins]
    %   Vf              - Array Transfer Functions [Q x nBins x K]
    %   Hf              - Target HRTFs [2 x nBins x K]
    %   f_sig           - Frequency vector [nBins x 1]
    %   lebedev_weights - Quadrature weights [K x 1] (Normalized sum to 1)
    
    
    Vf = permute(Vf,[1,3,2]);
    Hf = permute(Hf,[3,2,1]);

    [Q, nBins, K] = size(Vf);
    w = lebedev_weights(:) / sum(lebedev_weights); % Ensure normalized column vector
    
    % --- 1. Calculate Reference IC ---
    IC_ref = zeros(nBins, 1);
    for f = 1:nBins
        H_k = squeeze(Hf(:, f, :)); % [2 x K]
        C_ref = H_k * (w .* H_k');  % [2 x 2] Spatial Covariance
        
        c_ll = real(C_ref(1,1));
        c_rr = real(C_ref(2,2));
        c_lr = C_ref(1,2);
        
        IC_ref(f) = real(c_lr) / (sqrt(c_ll * c_rr) + 1e-12);
    end
    
    % --- 2. Calculate Estimated ICs ---
    IC_est = struct();
    for i = 1:length(algos_to_plot)
        name = algos_to_plot{i};
        if ~isfield(coeffs, name), continue; end
        
        W_algo = coeffs.(name).f_mls; % [Q x 2 x nBins]
        W_algo = permute(W_algo,[1,3,2]);
        IC_algo = zeros(nBins, 1);
        
        for f = 1:nBins
            W_f = squeeze(W_algo(:, :, f)); % [Q x 2]
            V_k = squeeze(Vf(:, f, :));     % [Q x K]
            
            % Synthesize the HRTF for all directions
            H_est = W_f' * V_k;             % [2 x K]
            
            % Compute Spatial Covariance
            C_est = H_est * (w .* H_est');  % [2 x 2]
            
            c_ll = real(C_est(1,1));
            c_rr = real(C_est(2,2));
            c_lr = C_est(1,2);
            
            IC_algo(f) = real(c_lr) / (sqrt(c_ll * c_rr) + 1e-12);
        end
        IC_est.(name) = IC_algo;
    end
    
    % --- 3. Plotting ---
    figure('Color', 'w', 'Name', 'Spatial Diffuse IC Analysis');
    plot(f_sig, IC_ref, 'k', 'LineWidth', 2, 'DisplayName', 'Reference (HRTF)');
    hold on;
    
    colors = lines(length(algos_to_plot));
    for i = 1:length(algos_to_plot)
        name = algos_to_plot{i};
        if isfield(IC_est, name)
            plot(f_sig, IC_est.(name), 'LineWidth', 1.5, 'Color', colors(i,:), 'DisplayName', name);
        end
    end
    
    set(gca, 'XScale', 'log');
    grid on;
    xlim([20, 20e3]);
    ylim([-0.2, 1.1]);
    xlabel('Frequency (Hz)');
    ylabel('Interaural Coherence (IC)');
    title('Spatial Diffuse IC (Integrated over Lebedev Grid)');
    legend('Location', 'southwest');
end

function c_BSM_constrained = apply_covariance_constraint_lebedev(c_BSM, ATF, HRTF, lebedev_weights)
    % Applies the Diffuse-Field Covariance Constraint using a Lebedev grid.
    % Inputs:
    %   c_BSM           - [Q x 2 x nBins] Original BSM weights
    %   ATF             - [Q x nBins x K] Array Transfer Functions 
    %   HRTF            - [2 x nBins x K] Head Related Transfer Functions
    %   lebedev_weights - [K x 1] Quadrature weights for the 1730 directions
    
    c_BSM = permute(c_BSM,[1,3,2]);
    ATF = permute(ATF,[1,3,2]);
    HRTF = permute(HRTF,[3,2,1]);

    [Q, nEars, nBins] = size(c_BSM);
    
    % 1. Format and Normalize the weights
    % Ensure it is a column vector
    w = lebedev_weights(:); 
    % Normalize so the integral over the sphere equals 1
    w = w / sum(w); 
    
    c_BSM_constrained = zeros(size(c_BSM));
    
    for f = 1:nBins
        % Extract data for this frequency bin
        W_f = squeeze(c_BSM(:, :, f)); % [Q x 2]
        V_f = squeeze(ATF(:, f, :));   % [Q x K]
        H_f = squeeze(HRTF(:, f, :));  % [2 x K]
        
        % 2. Calculate Target Diffuse Covariance (HRTF) using Quadrature
        % Fast matrix multiplication equivalent to SUM(w_k * h_k * h_k')
        C_target = H_f * (w .* H_f'); 
        
        % 3. Calculate Rendered Diffuse Covariance (BSM) using Quadrature
        R_v = V_f * (w .* V_f'); 
        C_rendered = W_f' * R_v * W_f;
        
        % Regularization to prevent singular matrix inversion
        beta = 1e-8 * trace(C_rendered);
        C_rendered_reg = C_rendered + beta * eye(2);
        
        % 4. Calculate Mixing Matrix M (Principal Matrix Square Roots)
        M = sqrtm(C_target) * inv(sqrtm(C_rendered_reg));
        
        % 5. Apply to Weights
        W_constrained = W_f * M';
        
        c_BSM_constrained(:, :, f) = W_constrained;
    end
    c_BSM_constrained = permute(c_BSM_constrained,[1,3,2]);
end


function c_BSM_constrained = apply_covariance_constraint(c_BSM, ATF, HRTF)
    % Applies the Diffuse-Field Covariance Constraint (Zotter & Frank, App A.5)
    % Inputs:
    %   c_BSM - [Q x 2 x nBins] Original BSM weights
    %   ATF   - [Q x nBins x K] Array Transfer Functions (Simulated/Measured)
    %   HRTF  - [2 x nBins x K] Head Related Transfer Functions
    % Output:
    %   c_BSM_constrained - [Q x 2 x nBins] Covariance-matched weights
    
    c_BSM = permute(c_BSM,[1,3,2]);
    ATF = permute(ATF,[1,3,2]);
    HRTF = permute(HRTF,[3,2,1]);
    [Q, nEars, nBins] = size(c_BSM);
    K = size(ATF, 3); % Number of directions in the diffuse field grid
    
    c_BSM_constrained = zeros(size(c_BSM));
    
    for f = 1:nBins
        % Extract data for this frequency bin
        W_f = squeeze(c_BSM(:, :, f)); % [Q x 2]
        V_f = squeeze(ATF(:, f, :));   % [Q x K]
        H_f = squeeze(HRTF(:, f, :));  % [2 x K]
        
        % 1. Calculate Target Diffuse Covariance (HRTF)
        % (Assuming grid is roughly uniform. Add quadrature weights if Lebedev)
        C_target = (H_f * H_f') / K; 
        
        % 2. Calculate Rendered Diffuse Covariance (BSM)
        R_v = (V_f * V_f') / K; % Spatial correlation of the array
        C_rendered = W_f' * R_v * W_f;
        
        % Regularization to prevent singular matrix inversion at low frequencies
        beta = 1e-8 * trace(C_rendered);
        C_rendered_reg = C_rendered + beta * eye(2);
        
        % 3. Calculate Mixing Matrix M (Matrix Square Roots)
        % Using sqrtm to find the principal square root
        M = sqrtm(C_target) * inv(sqrtm(C_rendered_reg));
        
        % 4. Apply to Weights
        % y_new = M * y_old = M * W^H * x = (W * M^H)^H * x
        W_constrained = W_f * M';
        
        c_BSM_constrained(:, :, f) = W_constrained;
    end
    c_BSM_constrained = permute(c_BSM_constrained,[1,3,2]);
end
    
function params = init_params(beta,arraysetup, base_path, save_folder_path)
%INIT_PARAMS Define all configurable parameters

    % Flags
    params.plot_figures = true;

    params.beta   = beta;
    params.base   = base_path;
    params.save   = base_path + save_folder_path + "/beta_" + num2str(beta) + "/";

    

    % Angles & array
    %params.angles       = [40,35]; % RBM POV (Vertical)
    params.angles       = [15,60]; % RBM POV (Vertical)
    params.alpha        = 1;
    %params.arrayType    = 1;
    %params.M            = 5;
    %params.arrayType    = 1;
    params.arrayType    = arraysetup.arrayType;
    params.M            = arraysetup.M;
    params.head_rot_az  = 0;

    % Source distribution
    params.source_distribution = 0;
    params.Q                   = 36;
    %params.f_cut_magLS         = 1.5e3;
    params.f_cut_magLS         = 3e3;

    % Paths
    params.HRTFpath    = base_path + "/HRTF/HRIR_L2702_koln.sofa";
    params.WignerDpath = base_path + "/third_party/WignerD/WignerDMatrix_diagN=32.mat";
    params.ATFpath     = base_path + "/ATF/aria_interpolated.mat";
end

function data = load_and_prepare_data(params)
%LOAD_AND_PREPARE_DATA Load HRTF/ATF, zero-pad, define grids

    startup_script(params.base)

    % Load HRTF
    hobj       = sofaToEaro(convertStringsToChars(params.HRTFpath));
    omega_hrtf = [hobj.sourceGrid.elevation,hobj.sourceGrid.azimuth,hobj.sourceGrid.r];

    % Load ATF
    ATF = load(params.ATFpath).aria;
    ATF.IR = ATF.IR_interpolated;
    ATF.IR = permute(ATF.IR,[1,3,2]);
    omega_atf = [ATF.lebedev_elevation,ATF.lebedev_azimuth]; %#ok<NASGU>

    % Azimuth plane grid
    az_res = 1;
    phi_az = deg2rad(linspace(-180,180,360/az_res).');

    %phi_az = deg2rad(linspace(0,180,180/az_res).');
    th_az  = deg2rad(90*ones(size(phi_az)));
    r_az   = ones(size(phi_az));
    omega_az = [th_az, phi_az, r_az];

    % Zero-pad HRTF and ATF to common nFFT
    tmp1 = 2^(nextpow2(size(ATF.IR,3)));
    tmp2 = 2^(nextpow2(size(hobj.data,2)));
    nFFT = max(tmp1,tmp2);
    %nFFT = 2^10;


    ATF.IR = cat(3, ATF.IR, zeros(size(ATF.IR,1),size(ATF.IR,2),nFFT-size(ATF.IR,3)));
    hobj.data = cat(2,hobj.data,zeros(size(hobj.data,1),nFFT-size(hobj.data,2),size(hobj.data,3)));
    hobj.taps = nFFT;

    % Store
    data.hobj = hobj;
    data.ATF = ATF;
    data.omega_hrtf = omega_hrtf;
    data.omega_az   = omega_az;
    data.nFFT = nFFT;
    data.fs   = hobj.fs;
end

function [objs,data] = setup_BSM_objects(params, data)
%SETUP_BSM_OBJECTS Create BSM objects with different distributions

    filt_len = data.nFFT / data.fs;
    M = params.M;
    at = params.arrayType;
    az = params.head_rot_az;
    fs = data.fs;

    objs.BSM    = setup_BSM_object(M,at,az,fs,filt_len,4,36,params.f_cut_magLS,params.angles);
    objs.FoV    = setup_BSM_object(M,at,az,fs,filt_len,3,36,params.f_cut_magLS,params.angles);
    objs.normal = setup_BSM_object(M,at,az,fs,filt_len,0,36,params.f_cut_magLS,params.angles);
    
    objs.BSM.beta = params.beta;

    data.omega_bsm        = [objs.BSM.th_BSMgrid_vec,objs.BSM.ph_BSMgrid_vec,ones(size(objs.BSM.ph_BSMgrid_vec))];
    data.omega_bsm_FoV    = [objs.FoV.th_BSMgrid_vec,objs.FoV.ph_BSMgrid_vec,ones(size(objs.FoV.ph_BSMgrid_vec))];
    data.omega_bsm_normal = [objs.normal.th_BSMgrid_vec,objs.normal.ph_BSMgrid_vec,ones(size(objs.normal.ph_BSMgrid_vec))];

    % Calculate sterring matrics
    [data.V_k_bsm, ~]        = get_me_some_V(objs.BSM,data.omega_bsm);
    [data.V_k_bsm_normal, ~] = get_me_some_V(objs.normal,data.omega_bsm_normal);
    [data.V_k_az, ~]         = get_me_some_V(objs.BSM,data.omega_az);

end

function [coeffs,data] = compute_BSM_filters(objs, data, params)
%COMPUTE_BSM_FILTERS Compute filter coefficients for each variant
    disp("COMPUTE_BSM_FILTERS Compute filter coefficients for each variant")

    W  = buildWeightedW([objs.BSM.th_BSMgrid_vec,objs.BSM.ph_BSMgrid_vec,ones(size(objs.BSM.ph_BSMgrid_vec))], ...
                        filterAzEl(data.omega_bsm, params.angles), params.alpha, params.beta);

    WJ = W; % (currently same)
    W_normal = eye(size([objs.normal.th_BSMgrid_vec,objs.normal.ph_BSMgrid_vec],1));

    data.W = W;
    data.WJ = WJ;
    data.W_normal = W_normal;
    % Pre process HRTF earo objects
    hobj_bsm        = process_hrtf(objs.BSM,params.HRTFpath,params.WignerDpath);
    hobj_bsm_normal = process_hrtf(objs.normal,params.HRTFpath,params.WignerDpath);


    % Compute coefficients
    % regular bsm
    disp("BSM")
    coeffs.BSM   = get_BSM_filters(objs.normal, hobj_bsm_normal, data.V_k_bsm_normal, W_normal);

    % tikoniv bsm
    disp("TBSM")
    objs.BSM.Tik = true;
    coeffs.TBSM  = get_BSM_filters(objs.BSM, hobj_bsm, data.V_k_bsm, W);
    objs.BSM.Tik = false;
    
    % wigted bsm
    disp("WBSM")
    coeffs.WBSM  = get_BSM_filters(objs.BSM, hobj_bsm, data.V_k_bsm, W);
    
    % Jan's bsm
    disp("JBSM")
    objs.BSM.Jan = true;
    coeffs.JBSM  = get_BSM_filters(objs.BSM, hobj_bsm, data.V_k_bsm, WJ);
    objs.BSM.Jan = false;


    % c_ls = coeffs.BSM.f_ls;
    % c_mls = coeffs.BSM.f_mls;
    % 
    % for f_idx = 1:data.nFFT/2 +1 
    %     c = c_ls(:,f_idx,1);
    %     V = data.V_k_az(:,:,f_idx);
    %     h_ls_hat(:,f_idx) = squeeze(c'*V);
    % end
    % h_ls_hat(:,end+1 : data.nFFT) = 0;
    % h_ls_hat_t = ifft(h_ls_hat, [], 2, 'symmetric');
    % 
    % for f_idx = 1:data.nFFT/2 +1 
    %     c = c_mls(:,f_idx,1);
    %     V = data.V_k_az(:,:,f_idx);
    %     h_mls_hat(:,f_idx) = squeeze(c'*V);
    % end
    % 
    % h_mls_hat(:,end+1 : data.nFFT) = 0;
    % h_mls_hat_t = ifft(h_mls_hat, [], 2, 'symmetric');
    % 
    % figure
    % plot_idx = 60;
    % plot(squeeze(h_ls_hat_t(plot_idx,:)))
    % hold on
    % plot(squeeze(h_mls_hat_t(plot_idx,:)))
    % legend(["LS","MagLS"])
    % c_ls_tmp = c_ls(:,:,1);
    % c_ls_tmp(:,end+1:data.nFFT,1)=0;
    % c_ls_tmp_time = ifft(c_ls_tmp, [], 2, 'symmetric');

end

function hrtf_hat = estimate_hrtf(objs, coeffs, data, params)
%ESTIMATE_HRTF Generate estimated HRTFs using filters (simplified)
    disp("ESTIMATE_HRTF Generate estimated HRTFs using filters")

    disp("calculating over a high lebedev grid")

    % data.V_k_bsm: [5 × 2702 × 257] = [coeff × dirs × freqs]
    % coeffs.*.f_* : [5 × 257 × 2]    = [coeff × freqs × ears]
    
    V = data.V_k_bsm;  % keep pages = freqs (3rd dim = 257)
    
    % Helper to avoid repeating permute:
    perm = @(X) permute(X, [1 3 2]);  % -> [5 × 2 × 257]
    
    % Each pagemtimes does: (coeffs^H [2×5×257]) * (V [5×2702×257]) = [2×2702×257]
    hrtf_hat.lebedev.ls.bsm  = pagemtimes(perm(coeffs.BSM.f_ls),  'ctranspose', V, 'none');
    hrtf_hat.lebedev.ls.wbsm = pagemtimes(perm(coeffs.WBSM.f_ls), 'ctranspose', V, 'none');
    hrtf_hat.lebedev.ls.tbsm = pagemtimes(perm(coeffs.TBSM.f_ls), 'ctranspose', V, 'none');
    hrtf_hat.lebedev.ls.jbsm = pagemtimes(perm(coeffs.JBSM.f_ls), 'ctranspose', V, 'none');
    
    hrtf_hat.lebedev.mls.bsm  = pagemtimes(perm(coeffs.BSM.f_mls),  'ctranspose', V, 'none');
    hrtf_hat.lebedev.mls.wbsm = pagemtimes(perm(coeffs.WBSM.f_mls), 'ctranspose', V, 'none');
    hrtf_hat.lebedev.mls.tbsm = pagemtimes(perm(coeffs.TBSM.f_mls), 'ctranspose', V, 'none');
    hrtf_hat.lebedev.mls.jbsm = pagemtimes(perm(coeffs.JBSM.f_mls), 'ctranspose', V, 'none');



    % disp("calculating over a high lebedev grid")
    % for idx_dir = 1:size(data.V_k_bsm,2)
    %     V_test = squeeze(data.V_k_bsm(:,idx_dir,:));
    %     for idx = 1:(data.nFFT/2+1)
    %         hrtf_hat.lebedev.ls.bsm(:,idx_dir,idx)  = squeeze(coeffs.BSM.f_ls(:,idx,:))'*V_test(:,idx);
    %         hrtf_hat.lebedev.ls.wbsm(:,idx_dir,idx) = squeeze(coeffs.WBSM.f_ls(:,idx,:))'*V_test(:,idx);
    %         hrtf_hat.lebedev.ls.tbsm(:,idx_dir,idx) = squeeze(coeffs.TBSM.f_ls(:,idx,:))'*V_test(:,idx);
    %         hrtf_hat.lebedev.ls.jbsm(:,idx_dir,idx) = squeeze(coeffs.JBSM.f_ls(:,idx,:))'*V_test(:,idx);
    % 
    %         hrtf_hat.lebedev.mls.bsm(:,idx_dir,idx)  = squeeze(coeffs.BSM.f_mls(:,idx,:))'*V_test(:,idx);
    %         hrtf_hat.lebedev.mls.wbsm(:,idx_dir,idx) = squeeze(coeffs.WBSM.f_mls(:,idx,:))'*V_test(:,idx);
    %         hrtf_hat.lebedev.mls.tbsm(:,idx_dir,idx) = squeeze(coeffs.TBSM.f_mls(:,idx,:))'*V_test(:,idx);
    %         hrtf_hat.lebedev.mls.jbsm(:,idx_dir,idx) = squeeze(coeffs.JBSM.f_mls(:,idx,:))'*V_test(:,idx);
    % 
    %     end
    % end

    disp("calculating over a dense horizontal grid")

    V = data.V_k_az;  % keep as [5 × nDirs × freqs]
    
    % Helper to permute coeffs into [5 × 2 × 257]
    perm = @(X) permute(X, [1 3 2]);
    
    % Each pagemtimes: (coeffs^H [2×5×257]) * (V [5×nDirs×257]) = [2×nDirs×257]
    
    hrtf_hat.az.ls.bsm  = pagemtimes(perm(coeffs.BSM.f_ls),  'ctranspose', V, 'none');
    hrtf_hat.az.ls.wbsm = pagemtimes(perm(coeffs.WBSM.f_ls), 'ctranspose', V, 'none');
    hrtf_hat.az.ls.tbsm = pagemtimes(perm(coeffs.TBSM.f_ls), 'ctranspose', V, 'none');
    hrtf_hat.az.ls.jbsm = pagemtimes(perm(coeffs.JBSM.f_ls), 'ctranspose', V, 'none');
    
    hrtf_hat.az.mls.bsm  = pagemtimes(perm(coeffs.BSM.f_mls),  'ctranspose', V, 'none');
    hrtf_hat.az.mls.wbsm = pagemtimes(perm(coeffs.WBSM.f_mls), 'ctranspose', V, 'none');
    hrtf_hat.az.mls.tbsm = pagemtimes(perm(coeffs.TBSM.f_mls), 'ctranspose', V, 'none');
    hrtf_hat.az.mls.jbsm = pagemtimes(perm(coeffs.JBSM.f_mls), 'ctranspose', V, 'none');


    % circshift the refference HRIR to the middle of the vector (lebedev)
    hobj_f_lebedev      = data.hobj;
    hobj_f_lebedev.data = circshift(hobj_f_lebedev.data,hobj_f_lebedev.taps/2 ,2);
    hobj_f_lebedev      = hobj_f_lebedev.toSH(44,'SRC');
    hobj_f_lebedev      = hobj_f_lebedev.toSpace('SRC',data.omega_bsm(:,1),data.omega_bsm(:,2));
    hobj_f_lebedev.data = real(hobj_f_lebedev.data);
    hobj_f_lebedev      = hobj_f_lebedev.toFreq(data.nFFT);
    hobj_f_lebedev.data = hobj_f_lebedev.data(:, 1:ceil(data.nFFT/2)+1, :);
    hrir_ref_lebedev    = ifft_pos_half(hobj_f_lebedev.data);

    % interpulate the horizontal HRIRs
    hobj_t_az = data.hobj;
    hobj_t_az.shutUp = true;
    hobj_t_az = hobj_t_az.toSH(44,'SRC');
    hobj_t_az = hobj_t_az.toSpace('SRC',data.omega_az(:,1),data.omega_az(:,2));
    hobj_t_az.data = real(hobj_t_az.data);
    hrir_ref_az = hobj_t_az.data;
    hrir_ref_az = circshift(hrir_ref_az,hobj_t_az.taps/2,2);


    disp("Aligning the estimations with the refference HRTF (lebedev)")
    is_freq = true;
    hrtf_hat.lebedev.ls.bsm   = post_process_estimates(hrtf_hat.lebedev.ls.bsm,hrir_ref_lebedev,is_freq);
    hrtf_hat.lebedev.ls.wbsm  = post_process_estimates(hrtf_hat.lebedev.ls.wbsm,hrir_ref_lebedev,is_freq);
    hrtf_hat.lebedev.ls.tbsm  = post_process_estimates(hrtf_hat.lebedev.ls.tbsm,hrir_ref_lebedev,is_freq);
    hrtf_hat.lebedev.ls.jbsm  = post_process_estimates(hrtf_hat.lebedev.ls.jbsm,hrir_ref_lebedev,is_freq);

    hrtf_hat.lebedev.mls.bsm   = post_process_estimates(hrtf_hat.lebedev.mls.bsm,hrir_ref_lebedev,is_freq);
    hrtf_hat.lebedev.mls.wbsm  = post_process_estimates(hrtf_hat.lebedev.mls.wbsm,hrir_ref_lebedev,is_freq);
    hrtf_hat.lebedev.mls.tbsm  = post_process_estimates(hrtf_hat.lebedev.mls.tbsm,hrir_ref_lebedev,is_freq);
    hrtf_hat.lebedev.mls.jbsm  = post_process_estimates(hrtf_hat.lebedev.mls.jbsm,hrir_ref_lebedev,is_freq);

    disp("Aligning the estimations with the refference HRTF (horizontal)")
    is_freq = false;
    hrtf_hat.az.ls.bsm    = post_process_estimates(hrtf_hat.az.ls.bsm,hrir_ref_az,is_freq);
    hrtf_hat.az.ls.wbsm   = post_process_estimates(hrtf_hat.az.ls.wbsm,hrir_ref_az,is_freq);
    hrtf_hat.az.ls.tbsm   = post_process_estimates(hrtf_hat.az.ls.tbsm,hrir_ref_az,is_freq);
    hrtf_hat.az.ls.jbsm   = post_process_estimates(hrtf_hat.az.ls.jbsm,hrir_ref_az,is_freq);

    hrtf_hat.az.mls.bsm    = post_process_estimates(hrtf_hat.az.mls.bsm,hrir_ref_az,is_freq);
    hrtf_hat.az.mls.wbsm   = post_process_estimates(hrtf_hat.az.mls.wbsm,hrir_ref_az,is_freq);
    hrtf_hat.az.mls.tbsm   = post_process_estimates(hrtf_hat.az.mls.tbsm,hrir_ref_az,is_freq);
    hrtf_hat.az.mls.jbsm   = post_process_estimates(hrtf_hat.az.mls.jbsm,hrir_ref_az,is_freq);


    hrtf_hat.az.ref = hrir_ref_az;
    % devide results to within FoV and outside FoV
    [inSector, ~] = ismember( data.omega_bsm, data.omega_bsm_FoV, 'rows');
    sector_idx    = find(inSector);
    outsector_idx = find(~inSector);

    hrtf_hat.lebedev.in.ref  = hobj_f_lebedev.data(sector_idx,:,:);
    hrtf_hat.lebedev.out.ref = hobj_f_lebedev.data(outsector_idx,:,:);

    hobj_f_lebedev_jan = zeros(size(hobj_f_lebedev.data));
    for ear_idx = 1:2
        %hobj_f_lebedev_jan(:,:,ear_idx) = data.WJ*hobj_f_lebedev.data(:,:,ear_idx);
        hobj_f_lebedev_jan(:,:,ear_idx) = hobj_f_lebedev.data(:,:,ear_idx);

    end
    hrtf_hat.lebedev.in.ref_jan  = hobj_f_lebedev_jan(sector_idx,:,:);
    hrtf_hat.lebedev.out.ref_jan = hobj_f_lebedev_jan(outsector_idx,:,:);

    hrtf_hat.lebedev.in.bsm.ls   = hrtf_hat.lebedev.ls.bsm(sector_idx,:,:);
    hrtf_hat.lebedev.out.bsm.ls  = hrtf_hat.lebedev.ls.bsm(outsector_idx,:,:);
    hrtf_hat.lebedev.in.bsm.mls  = hrtf_hat.lebedev.mls.bsm(sector_idx,:,:);
    hrtf_hat.lebedev.out.bsm.mls = hrtf_hat.lebedev.mls.bsm(outsector_idx,:,:);

    hrtf_hat.lebedev.in.wbsm.ls   = hrtf_hat.lebedev.ls.wbsm(sector_idx,:,:);
    hrtf_hat.lebedev.out.wbsm.ls  = hrtf_hat.lebedev.ls.wbsm(outsector_idx,:,:);
    hrtf_hat.lebedev.in.wbsm.mls  = hrtf_hat.lebedev.mls.wbsm(sector_idx,:,:);
    hrtf_hat.lebedev.out.wbsm.mls = hrtf_hat.lebedev.mls.wbsm(outsector_idx,:,:);

    hrtf_hat.lebedev.in.tbsm.ls   = hrtf_hat.lebedev.ls.tbsm(sector_idx,:,:);
    hrtf_hat.lebedev.out.tbsm.ls  = hrtf_hat.lebedev.ls.tbsm(outsector_idx,:,:);
    hrtf_hat.lebedev.in.tbsm.mls  = hrtf_hat.lebedev.mls.tbsm(sector_idx,:,:);
    hrtf_hat.lebedev.out.tbsm.mls = hrtf_hat.lebedev.mls.tbsm(outsector_idx,:,:);

    hrtf_hat.lebedev.in.jbsm.ls   = hrtf_hat.lebedev.ls.jbsm(sector_idx,:,:);
    hrtf_hat.lebedev.out.jbsm.ls  = hrtf_hat.lebedev.ls.jbsm(outsector_idx,:,:);
    hrtf_hat.lebedev.in.jbsm.mls  = hrtf_hat.lebedev.mls.jbsm(sector_idx,:,:);
    hrtf_hat.lebedev.out.jbsm.mls = hrtf_hat.lebedev.mls.jbsm(outsector_idx,:,:);

    hrtf_hat.lebedev = rmfield(hrtf_hat.lebedev, 'ls');
    hrtf_hat.lebedev = rmfield(hrtf_hat.lebedev, 'mls');

    hrtf_hat.az.bsm.ls  = hrtf_hat.az.ls.bsm;
    hrtf_hat.az.wbsm.ls = hrtf_hat.az.ls.wbsm;
    hrtf_hat.az.tbsm.ls = hrtf_hat.az.ls.tbsm;
    hrtf_hat.az.jbsm.ls = hrtf_hat.az.ls.jbsm;
    
    hrtf_hat.az.bsm.mls = hrtf_hat.az.mls.bsm;
    hrtf_hat.az.wbsm.mls = hrtf_hat.az.mls.wbsm;
    hrtf_hat.az.tbsm.mls = hrtf_hat.az.mls.tbsm;
    hrtf_hat.az.jbsm.mls = hrtf_hat.az.mls.jbsm;

    hrtf_hat.az = rmfield(hrtf_hat.az, 'ls');
    hrtf_hat.az = rmfield(hrtf_hat.az, 'mls');

    hrtf_hat.H_t  = hrir_ref_lebedev;
    hrtf_hat.H_f  = hobj_f_lebedev.data;
    
    group_delays_ref = group_delay_multi_channel(hrtf_hat.az.ref(:,:,1).',size(hrtf_hat.az.ref,2));
    group_delays_bsm_mls = group_delay_multi_channel(hrtf_hat.az.bsm.mls(:,:,1).',size(hrtf_hat.az.bsm.mls,2));
    group_delays_bsm_ls = group_delay_multi_channel(hrtf_hat.az.bsm.ls(:,:,1).',size(hrtf_hat.az.bsm.ls,2));

    % figure
    % subplot(3,1,1)
    % semilogx(group_delays_ref)
    % xlim([0 300])
    % ylim([500 540])
    % subplot(3,1,2)
    % semilogx(group_delays_bsm_mls)
    % xlim([0 300])
    % ylim([500 540])
    % subplot(3,1,3)
    % semilogx(group_delays_bsm_ls)
    % xlim([0 300])
    % ylim([500 540])


end

function hrir = ifft_pos_half(hrtf)
    % Calc ILD errors
    nFFT = 2*(size(hrtf,2)-1);
    H_ref_full = zeros(size(hrtf,1), nFFT, size(hrtf,3));
    H_ref_full(:, 1:nFFT/2+1, :) = hrtf;
    H_ref_full(:, nFFT/2+2:end, :) = conj(flip(hrtf(:, 2:nFFT/2, :), 2));
    hrir = ifft(H_ref_full, nFFT, 2, 'symmetric');

end

function h_out = post_process_estimates(h_in,h_ref,is_freq)
    target_lufs = -23;
    h_in = permute(h_in,[2,3,1]);
    nFFT = 2*(size(h_in,2)-1);
    
    % aligne in time    
    h_in_t = ifft_pos_half(h_in);
    
    %h_in_t = AlignFilters_Time(h_ref,h_in_t);
    if is_freq
        h_in_f = fft(h_in_t, nFFT, 2);
        h_out = h_in_f(:,1:nFFT/2+1,:);
    else
        h_out = h_in_t;
    end




end

function [errors,data] = evaluate_errors(hrtf, data, objs, params)
%EVALUATE_ERRORS Compute ILD/ITD/NMSE/NMAGE/BSD metrics (simplified)

    disp("EVALUATE_ERRORS Compute ILD/ITD/NMSE/NMAGE/BSD metrics")
    disp("Compute ILD")
    ild_f_band = [1.5e3, 8e3];
    [errors.ild.ls.REF,~]   = calc_ILD(hrtf.az.ref, ild_f_band, data.fs);
    errors.ild.mls.REF      = errors.ild.ls.REF;
    [errors.ild.ls.BSM,~]   = calc_ILD(hrtf.az.bsm.ls, ild_f_band, data.fs);
    [errors.ild.mls.BSM,~]  = calc_ILD(hrtf.az.bsm.mls, ild_f_band, data.fs);
    %[errors.ild.ls.WBSM,~]  = calc_ILD(hrtf.az.wbsm.ls, ild_f_band, data.fs);
    %[errors.ild.mls.WBSM,~] = calc_ILD(hrtf.az.wbsm.mls, ild_f_band, data.fs);
    [errors.ild.ls.WBSM,~]  = calc_ILD(hrtf.az.tbsm.ls, ild_f_band, data.fs);
    [errors.ild.mls.WBSM,data.f_c_ild] = calc_ILD(hrtf.az.tbsm.mls, ild_f_band, data.fs);
    %[errors.ild.ls.JBSM,~]  = calc_ILD(hrtf.az.jbsm.ls, ild_f_band, data.fs);
    %[errors.ild.mls.JBSM,data.f_c_ild] = calc_ILD(hrtf.az.jbsm.mls, ild_f_band, data.fs);

    disp("Compute ILD - on the sphere")
    [errors.ild_s.in.ls.REF,~]   = calc_ILD(hrtf.lebedev.in.ref, ild_f_band, data.fs);
    errors.ild_s.in.mls.REF      = errors.ild_s.in.ls.REF;
    [errors.ild_s.in.ls.BSM,~]   = calc_ILD(hrtf.lebedev.in.bsm.ls, ild_f_band, data.fs);
    [errors.ild_s.in.mls.BSM,~]  = calc_ILD(hrtf.lebedev.in.bsm.mls, ild_f_band, data.fs);
    %[errors.ild_s.in.ls.WBSM,~]  = calc_ILD(hrtf.lebedev.in.wbsm.ls, ild_f_band, data.fs);
    %[errors.ild_s.in.mls.WBSM,~] = calc_ILD(hrtf.lebedev.in.wbsm.mls, ild_f_band, data.fs);
    [errors.ild_s.in.ls.WBSM,~]  = calc_ILD(hrtf.lebedev.in.tbsm.ls, ild_f_band, data.fs);
    [errors.ild_s.in.mls.WBSM,~] = calc_ILD(hrtf.lebedev.in.tbsm.mls, ild_f_band, data.fs);
    %[errors.ild_s.in.ls.JBSM,~]  = calc_ILD(hrtf.lebedev.in.jbsm.ls, ild_f_band, data.fs);
    %[errors.ild_s.in.mls.JBSM,~] = calc_ILD(hrtf.lebedev.in.jbsm.mls, ild_f_band, data.fs);

    [errors.ild_s.out.ls.REF,~]   = calc_ILD(hrtf.lebedev.out.ref, ild_f_band, data.fs);
    errors.ild_s.out.mls.REF      = errors.ild_s.out.ls.REF;
    [errors.ild_s.out.ls.BSM,~]   = calc_ILD(hrtf.lebedev.out.bsm.ls, ild_f_band, data.fs);
    [errors.ild_s.out.mls.BSM,~]  = calc_ILD(hrtf.lebedev.out.bsm.mls, ild_f_band, data.fs);
    %[errors.ild_s.out.ls.WBSM,~]  = calc_ILD(hrtf.lebedev.out.wbsm.ls, ild_f_band, data.fs);
    %[errors.ild_s.out.mls.WBSM,~] = calc_ILD(hrtf.lebedev.out.wbsm.mls, ild_f_band, data.fs);
    [errors.ild_s.out.ls.WBSM,~]  = calc_ILD(hrtf.lebedev.out.tbsm.ls, ild_f_band, data.fs);
    [errors.ild_s.out.mls.WBSM,~] = calc_ILD(hrtf.lebedev.out.tbsm.mls, ild_f_band, data.fs);
    %[errors.ild_s.out.ls.JBSM,~]  = calc_ILD(hrtf.lebedev.out.jbsm.ls, ild_f_band, data.fs);
    %[errors.ild_s.out.mls.JBSM,~] = calc_ILD(hrtf.lebedev.in.jbsm.mls, ild_f_band, data.fs);



    disp("Compute ITD")
    %gd_cut     = 1.45e3;

    gd_cut     = 1.2e3;
    errors.itd.ls.REF   = calc_ITD(hrtf.az.ref, gd_cut, data.fs, objs.BSM.freqs_sig);
    errors.itd.mls.REF  = errors.itd.ls.REF;
    errors.itd.ls.BSM   = calc_ITD(hrtf.az.bsm.ls, gd_cut, data.fs, objs.BSM.freqs_sig);
    errors.itd.mls.BSM  = calc_ITD(hrtf.az.bsm.mls, gd_cut, data.fs, objs.BSM.freqs_sig);
    %errors.itd.ls.WBSM  = calc_ITD(hrtf.az.wbsm.ls, gd_cut, data.fs, objs.BSM.freqs_sig);
    %errors.itd.mls.WBSM = calc_ITD(hrtf.az.wbsm.mls, gd_cut, data.fs, objs.BSM.freqs_sig);
    errors.itd.ls.WBSM  = calc_ITD(hrtf.az.tbsm.ls, gd_cut, data.fs, objs.BSM.freqs_sig);
    errors.itd.mls.WBSM = calc_ITD(hrtf.az.tbsm.mls, gd_cut, data.fs, objs.BSM.freqs_sig);
    %errors.itd.ls.JBSM  = calc_ITD(hrtf.az.jbsm.ls, gd_cut, data.fs, objs.BSM.freqs_sig);
    %errors.itd.mls.JBSM = calc_ITD(hrtf.az.jbsm.mls, gd_cut, data.fs, objs.BSM.freqs_sig);

    disp("Compute NMSE")
    [errors.nmse.in.ls.BSM.mu,errors.nmse.in.ls.BSM.std]    = calc_nmse(hrtf.lebedev.in.bsm.ls,hrtf.lebedev.in.ref);
    [errors.nmse.in.mls.BSM.mu,errors.nmse.in.mls.BSM.std]  = calc_nmse(hrtf.lebedev.in.bsm.mls,hrtf.lebedev.in.ref);
    [errors.nmse.out.ls.BSM.mu,errors.nmse.out.ls.BSM.std]   = calc_nmse(hrtf.lebedev.out.bsm.ls,hrtf.lebedev.out.ref);
    [errors.nmse.out.mls.BSM.mu,errors.nmse.out.mls.BSM.std] = calc_nmse(hrtf.lebedev.out.bsm.mls,hrtf.lebedev.out.ref);

    %[errors.nmse.in.ls.WBSM.mu,errors.nmse.in.ls.WBSM.std]    = calc_nmse(hrtf.lebedev.in.wbsm.ls,hrtf.lebedev.in.ref);
    %[errors.nmse.in.mls.WBSM.mu,errors.nmse.in.mls.WBSM.std]  = calc_nmse(hrtf.lebedev.in.wbsm.mls,hrtf.lebedev.in.ref);
    %[errors.nmse.out.ls.WBSM.mu,errors.nmse.out.ls.WBSM.std]   = calc_nmse(hrtf.lebedev.out.wbsm.ls,hrtf.lebedev.out.ref);
    %[errors.nmse.out.mls.WBSM.mu,errors.nmse.out.mls.WBSM.std] = calc_nmse(hrtf.lebedev.out.wbsm.mls,hrtf.lebedev.out.ref);

    [errors.nmse.in.ls.WBSM.mu,errors.nmse.in.ls.WBSM.std]    = calc_nmse(hrtf.lebedev.in.tbsm.ls,hrtf.lebedev.in.ref);
    [errors.nmse.in.mls.WBSM.mu,errors.nmse.in.mls.WBSM.std]  = calc_nmse(hrtf.lebedev.in.tbsm.mls,hrtf.lebedev.in.ref);
    [errors.nmse.out.ls.WBSM.mu,errors.nmse.out.ls.WBSM.std]   = calc_nmse(hrtf.lebedev.out.tbsm.ls,hrtf.lebedev.out.ref);
    [errors.nmse.out.mls.WBSM.mu,errors.nmse.out.mls.WBSM.std] = calc_nmse(hrtf.lebedev.out.tbsm.mls,hrtf.lebedev.out.ref);

    %[errors.nmse.in.ls.JBSM.mu,errors.nmse.in.ls.JBSM.std]    = calc_nmse(hrtf.lebedev.in.jbsm.ls,hrtf.lebedev.in.ref);
    %[errors.nmse.in.mls.JBSM.mu,errors.nmse.in.mls.JBSM.std]  = calc_nmse(hrtf.lebedev.in.jbsm.mls,hrtf.lebedev.in.ref);
    %[errors.nmse.out.ls.JBSM.mu,errors.nmse.out.ls.JBSM.std]   = calc_nmse(hrtf.lebedev.out.jbsm.ls,hrtf.lebedev.out.ref);
    %[errors.nmse.out.mls.JBSM.mu,errors.nmse.out.mls.JBSM.std] = calc_nmse(hrtf.lebedev.out.jbsm.mls,hrtf.lebedev.out.ref);

    disp("Compute Normalized Magntidue Error")
    [errors.mag.in.ls.BSM.mu,errors.mag.in.ls.BSM.std]    = calc_mag(hrtf.lebedev.in.bsm.ls,hrtf.lebedev.in.ref);
    [errors.mag.in.mls.BSM.mu,errors.mag.in.mls.BSM.std]  = calc_mag(hrtf.lebedev.in.bsm.mls,hrtf.lebedev.in.ref);
    [errors.mag.out.ls.BSM.mu,errors.mag.out.ls.BSM.std]   = calc_mag(hrtf.lebedev.out.bsm.ls,hrtf.lebedev.out.ref);
    [errors.mag.out.mls.BSM.mu,errors.mag.out.mls.BSM.std] = calc_mag(hrtf.lebedev.out.bsm.mls,hrtf.lebedev.out.ref);

    %[errors.mag.in.ls.WBSM.mu,errors.mag.in.ls.WBSM.std]    = calc_mag(hrtf.lebedev.in.wbsm.ls,hrtf.lebedev.in.ref);
    %[errors.mag.in.mls.WBSM.mu,errors.mag.in.mls.WBSM.std]  = calc_mag(hrtf.lebedev.in.wbsm.mls,hrtf.lebedev.in.ref);
    %[errors.mag.out.ls.WBSM.mu,errors.mag.out.ls.WBSM.std]   = calc_mag(hrtf.lebedev.out.wbsm.ls,hrtf.lebedev.out.ref);
    %[errors.mag.out.mls.WBSM.mu,errors.mag.out.mls.WBSM.std] = calc_mag(hrtf.lebedev.out.wbsm.mls,hrtf.lebedev.out.ref);

    [errors.mag.in.ls.WBSM.mu,errors.mag.in.ls.WBSM.std]    = calc_mag(hrtf.lebedev.in.tbsm.ls,hrtf.lebedev.in.ref);
    [errors.mag.in.mls.WBSM.mu,errors.mag.in.mls.WBSM.std]  = calc_mag(hrtf.lebedev.in.tbsm.mls,hrtf.lebedev.in.ref);
    [errors.mag.out.ls.WBSM.mu,errors.mag.out.ls.WBSM.std]   = calc_mag(hrtf.lebedev.out.tbsm.ls,hrtf.lebedev.out.ref);
    [errors.mag.out.mls.WBSM.mu,errors.mag.out.mls.WBSM.std] = calc_mag(hrtf.lebedev.out.tbsm.mls,hrtf.lebedev.out.ref);

    %[errors.mag.in.ls.JBSM.mu,errors.mag.in.ls.JBSM.std]    = calc_mag(hrtf.lebedev.in.jbsm.ls,hrtf.lebedev.in.ref);
    %[errors.mag.in.mls.JBSM.mu,errors.mag.in.mls.JBSM.std]  = calc_mag(hrtf.lebedev.in.jbsm.mls,hrtf.lebedev.in.ref);
    %[errors.mag.out.ls.JBSM.mu,errors.mag.out.ls.JBSM.std]   = calc_mag(hrtf.lebedev.out.jbsm.ls,hrtf.lebedev.out.ref);
    %[errors.mag.out.mls.JBSM.mu,errors.mag.out.mls.JBSM.std] = calc_mag(hrtf.lebedev.out.jbsm.mls,hrtf.lebedev.out.ref);

    disp("Compute BSD errors")
    f_band = [20 , data.fs/2];
    ERBstruct = prepare_ERB(f_band, data.fs);


    [errors.bsd.in.ls.BSM, data.f_c_bsd]   = get_gammatone(hrtf.lebedev.in.ref,hrtf.lebedev.in.bsm.ls,f_band,data.fs,ERBstruct);
    [errors.bsd.in.mls.BSM,~]  = get_gammatone(hrtf.lebedev.in.ref,hrtf.lebedev.in.bsm.mls,f_band,data.fs,ERBstruct);
    [errors.bsd.out.ls.BSM,~]  = get_gammatone(hrtf.lebedev.out.ref,hrtf.lebedev.out.bsm.ls,f_band,data.fs,ERBstruct);
    [errors.bsd.out.mls.BSM,~] = get_gammatone(hrtf.lebedev.out.ref,hrtf.lebedev.out.bsm.mls,f_band,data.fs,ERBstruct);

    %[errors.bsd.in.ls.WBSM,~]   = get_gammatone(hrtf.lebedev.in.ref,hrtf.lebedev.in.wbsm.ls,f_band,data.fs,ERBstruct);
    %[errors.bsd.in.mls.WBSM,~]  = get_gammatone(hrtf.lebedev.in.ref,hrtf.lebedev.in.wbsm.mls,f_band,data.fs,ERBstruct);
    %[errors.bsd.out.ls.WBSM,~]  = get_gammatone(hrtf.lebedev.out.ref,hrtf.lebedev.out.wbsm.ls,f_band,data.fs,ERBstruct);
    %[errors.bsd.out.mls.WBSM,~] = get_gammatone(hrtf.lebedev.out.ref,hrtf.lebedev.out.wbsm.mls,f_band,data.fs,ERBstruct);

    [errors.bsd.in.ls.WBSM,~]   = get_gammatone(hrtf.lebedev.in.ref,hrtf.lebedev.in.tbsm.ls,f_band,data.fs,ERBstruct);
    [errors.bsd.in.mls.WBSM,~]  = get_gammatone(hrtf.lebedev.in.ref,hrtf.lebedev.in.tbsm.mls,f_band,data.fs,ERBstruct);
    [errors.bsd.out.ls.WBSM,~]  = get_gammatone(hrtf.lebedev.out.ref,hrtf.lebedev.out.tbsm.ls,f_band,data.fs,ERBstruct);
    [errors.bsd.out.mls.WBSM,~] = get_gammatone(hrtf.lebedev.out.ref,hrtf.lebedev.out.tbsm.mls,f_band,data.fs,ERBstruct);

    %[errors.bsd.in.ls.JBSM,~]   = get_gammatone(hrtf.lebedev.in.ref,hrtf.lebedev.in.jbsm.ls,f_band,data.fs,ERBstruct);
    %[errors.bsd.in.mls.JBSM,~]  = get_gammatone(hrtf.lebedev.in.ref,hrtf.lebedev.in.jbsm.mls,f_band,data.fs,ERBstruct);
    %[errors.bsd.out.ls.JBSM,~]  = get_gammatone(hrtf.lebedev.out.ref,hrtf.lebedev.out.jbsm.ls,f_band,data.fs,ERBstruct);
    %[errors.bsd.out.mls.JBSM,~] = get_gammatone(hrtf.lebedev.out.ref,hrtf.lebedev.out.jbsm.mls,f_band,data.fs,ERBstruct);


    errors = evaluate_all_dirs_errors(errors, hrtf);

end

function [nmse, nmse_std] = calc_mag(x, y)

    % Compute squared error
    error_squared = abs(abs(x) - abs(y)).^2;

    % Compute mean square error across directions
    mse_dir = mean(error_squared, 1);  % size: [1 x freq_bins x ears]
    std_dir = std(error_squared, 0, 1); % std across directions

    % Compute power of ground truth (normalization factor)
    power_ref = mean(abs(y).^2, 1);  % same shape
    power_std = std(abs(y).^2, 0, 1); % optional, not used

    % Avoid divide-by-zero
    eps_val = 1e-12;
    nmse = mse_dir ./ (power_ref + eps_val);  % normalized MSE
    nmse_std = std_dir ./ (power_ref + eps_val);  % std of NMSE

    % Squeeze to get [freq_bins x ears]
    nmse = squeeze(nmse);
    nmse_std = squeeze(nmse_std);
end




function [nmse, nmse_std] = calc_nmse(x, y)

    % Compute squared error
    error_squared = abs(x - y).^2;

    % Compute mean square error across directions
    mse_dir = mean(error_squared, 1);  % size: [1 x freq_bins x ears]
    std_dir = std(error_squared, 0, 1); % std across directions

    % Compute power of ground truth (normalization factor)
    power_ref = mean(abs(y).^2, 1);  % same shape
    power_std = std(abs(y).^2, 0, 1); % optional, not used

    % Avoid divide-by-zero
    eps_val = 1e-12;
    nmse = mse_dir ./ (power_ref + eps_val);  % normalized MSE
    nmse_std = std_dir ./ (power_ref + eps_val);  % std of NMSE

    % Squeeze to get [freq_bins x ears]
    nmse = squeeze(nmse);
    nmse_std = squeeze(nmse_std);
end


function alpha = get_optimal_gain(x, y)
    % Find the optimal real scalar 'alpha' to minimize MSE(y, alpha * x)
    % x: estimated signal (HRTF_hat or HRIR_hat)
    % y: reference signal (HRTF or HRIR)
    
    % Flatten the arrays to 1D vectors to sum across all dimensions
    x_flat = x(:);
    y_flat = y(:);
    
    % Numerator: Sum of the real part of the inner product
    numerator = sum(real(conj(x_flat) .* y_flat));
    
    % Denominator: Total energy of the estimate
    denominator = sum(abs(x_flat).^2);
    
    % Avoid division by zero in case the estimate is pure silence
    eps_val = 1e-12; 
    
    % Calculate the global gain factor
    alpha = numerator / (denominator + eps_val);
end

function [ILD_no_av,f_c] = calc_ILD(P_t,f_band,fs)
    P_t_l = P_t(:,:,1).';
    P_t_r = P_t(:,:,2).';
    
    [ILD_no_av,f_c,~] = AKerbILD(P_t_l, P_t_r, f_band, fs);
end

function ITD = calc_ITD(H_t,gd_cut,fs,f_vec)

    [~,gd_cut_idx] = min(abs(f_vec- gd_cut));
    nfft = (length(f_vec)-1)*2;

    
    group_delays_L = group_delay_multi_channel(H_t(:,:,1).',nfft);
    group_delays_R = group_delay_multi_channel(H_t(:,:,2).',nfft);
    
    TOA_L = mean(group_delays_L(1:gd_cut_idx,:),1);
    TOA_R = mean(group_delays_R(1:gd_cut_idx,:),1);
    ITD = TOA_R - TOA_L; % in samples
    ITD = (ITD/fs)*1e6;  % in microsecs
    
end

function group_delays = group_delay_multi_channel(impulse_response, N_fft)


    % group_delay_multi_channel Compute group delay for multi-channel impulse responses
    %
    % Inputs:
    %   impulse_response : [N x ch] time-domain impulse responses
    %   N_fft (optional) : number of FFT points (zero-padding if N_fft > N)
    %
    % Output:
    %   group_delays     : [N_fft/2+1 x ch] group delay per frequency bin and channel

    if nargin < 2
        N_fft = size(impulse_response, 1);  % No padding by default
    end

    [N, ch] = size(impulse_response);

    % Time index vector
    n = (0:N-1)';
    nh = impulse_response .* n;

    % Compute FFT with zero-padding to N_fft
    dft_h = fft(impulse_response, N_fft, 1) / N;
    dft_nh = fft(nh, N_fft, 1) / N;

    % Compute group delay
    eps_val = 1e-8;
    group_delays = real(dft_nh ./ (dft_h + eps_val));

    % Keep only the positive frequency bins
    group_delays = group_delays(1:floor(N_fft/2)+1, :);

    % Check for NaNs
    if any(isnan(group_delays), 'all')
        warning('NaNs detected in group delay!');
        min_val = min(abs(dft_h), [], 'all');
        fprintf('Min abs(dft_h): %g\n', min_val);
        error('NaNs in group delay computation');
    end
end

function plot_results(errors, hrtf,data ,objs, params)
%PLOT_RESULTS Generate plots
    title_name = "Binaural Erros";
    plot_combined_errors(objs.BSM, errors,data.f_c_bsd, title_name, params.save);

    algorithms_to_plot = {'BSM', 'TBSM'}; 
    plot_all_dirs_errors(objs.BSM, errors, title_name, params.save, algorithms_to_plot);

    e_ls = errors.ild.ls;
    e_mls = errors.ild.mls;
    plot_ILD_curves(e_ls, rad2deg(data.omega_az(:,2)), "LS ILD curves", params.save);
    plot_ILD_curves(e_mls, rad2deg(data.omega_az(:,2)), "MagLS ILD curves", params.save);
    

    plot_ILD_error_curves(e_ls, rad2deg(data.omega_az(:,2)),data.f_c_ild, "LS ", params.save);
    plot_ILD_error_curves(e_mls, rad2deg(data.omega_az(:,2)),data.f_c_ild, "MagLS ", params.save);

    e_ls = errors.itd.ls;
    e_mls = errors.itd.mls;
    plot_ITD_curves(e_ls, rad2deg(data.omega_az(:,2)), "LS ITD curves", params.save);
    plot_ITD_curves(e_mls, rad2deg(data.omega_az(:,2)), "MagLS ITD curves", params.save);

    % title_name = "freq_";
    % plot_errors_freq(objs.BSM,errors,title_name,params.save)
    % 
    % title_name = "ERB-BSD";
    % plot_errors_freq_BSD(errors,data.f_c_bsd,title_name,params.save);



end

function plot_combined_errors(BSMobj, e, fc, title_name, save_folder)
% PLOT_COMBINED_ERRORS Generates a 2x2 grid for conference papers:
% Top row: NMSE (In/Out), Bottom row: BSD (In/Out)
%
% Inputs:
%   - BSMobj: Object with .freqs_sig (for NMSE)
%   - e: Struct containing .nmse and .bsd nested fields
%   - fc: Frequency axis for BSD (ERB scale)
%   - title_name: Base string for saving
%   - save_folder: Destination path

    algos = fieldnames(e.nmse.in.ls);
    colors = lines(numel(algos));
    f_sig = BSMobj.freqs_sig.'; % Frequency axis for NMSE
    f_erb = fc(:);              % Frequency axis for BSD (ensured as column)
    
    optimization = {'ls', 'mls'}; 

    for opt_idx = 1:numel(optimization)
        opt_type = optimization{opt_idx};
        
        fig = figure; 
        t = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
        
        % Plot order: [NMSE_In, NMSE_Out, BSD_In, BSD_Out]
        metrics = {'nmse', 'nmse', 'mag', 'mag', 'bsd', 'bsd'};
        dirs = {'in', 'out', 'in', 'out', 'in', 'out'};
        titles = {'NMSE - Within FoV', 'NMSE - Outside FoV', ...
                  'MagE - Within FoV', 'MagE - Outside FoV', ...
                  'BSD - Within FoV', 'BSD - Outside FoV'};

        for i = 1:length(metrics)
            ax = nexttile;
            hold on;
            curr_metric = metrics{i};
            curr_dir = dirs{i};
            
            % Select the correct frequency vector for the current metric
            if strcmp(curr_metric, 'bsd')
                f_curr = f_erb;
            else
                f_curr = f_sig;
            end
            
            for k = 1:numel(algos)
                name = algos{k};
                
                if strcmp(curr_metric, 'nmse')
                    % NMSE Logic
                    mu_lin  = mean(e.nmse.(curr_dir).(opt_type).(name).mu, 2);
                    std_lin = mean(e.nmse.(curr_dir).(opt_type).(name).std, 2);
                    
                    mu_plot = 10*log10(mu_lin);
                    sigma   = 10*log10(mu_lin + std_lin) - mu_plot;
                    y_label_str = 'NMSE (dB)';
                    y_lims = [-40, 10];
                elseif strcmp(curr_metric, 'mag')
                    % Normalized Magnitude Logic
                    mu_lin  = mean(e.mag.(curr_dir).(opt_type).(name).mu, 2);
                    std_lin = mean(e.mag.(curr_dir).(opt_type).(name).std, 2);
                    
                    mu_plot = 10*log10(mu_lin);
                    sigma   = 10*log10(mu_lin + std_lin) - mu_plot;
                    y_label_str = 'MagError (dB)';
                    y_lims = [-40, 10];
                    
                    
                else
                    % BSD Logic (Averaging across ears)
                    bsd_raw = mean(e.bsd.(curr_dir).(opt_type).(name), 3); 
                    mu_plot = mean(bsd_raw, 2);
                    sigma   = std(bsd_raw, 0, 2);
                    y_label_str = 'BSD (dB)';
                    y_lims = [0, 7.5];
                end

                % Shaded Error Region
                % Using f_curr to ensure lengths match the data (mu_plot/sigma)
                fill([f_curr; flipud(f_curr)], [mu_plot - sigma; flipud(mu_plot + sigma)], ...
                     colors(k,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

                % Mean Line
                plot(f_curr, mu_plot, 'LineWidth', 1.2, 'Color', colors(k,:));
            end

            % Formatting
            set(ax, 'XScale', 'log');
            grid on;
            xlim([20, 20e3]);
            ylim(y_lims);
            title(titles{i});
            
            if mod(i, 3) ~= 0, ylabel(y_label_str); end
            if i > 3, xlabel('Frequency (Hz)'); end
        end

        % Global Legend
        lgd = legend(algos, 'Orientation', 'horizontal', 'Location', 'southoutside');
        lgd.Layout.Tile = 'south';

        % Save
        fig_name = fullfile(save_folder, sprintf('%s_%s_combined.fig', title_name, opt_type));
        savefig(fig, fig_name);
    end
end

function plot_errors_freq_BSD(e, fc, title_name, save_folder)
%PLOT_ERRORS_FREQ_BSD Plots BSD error curves (single figure with all algorithms)
%   Each of the 4 algorithms appears with two lines: solid for Within FoV and dashed for Outside FoV
%
%   Inputs:
%     - e: struct with fields e.bsd.{BSM,WBSM,TBSM,CBSM}.{in,out}.mu (size [nFFT x 2])
%     - fc: frequency axis [1 x nFFT]
%     - title_name: plot title (string)
%     - save_folder: folder to save .fig file (string)

    algos = fieldnames(e.bsd.in.ls);
    colors = lines(numel(algos));
    dirs     = {'in','out'};
    optimization = {'ls','mls'};
    dirName  = {'Within','Outside'};


    for opt_curr = 1:numel(optimization)
        opt_curr_n = optimization{opt_curr};
        for d = 1:numel(dirs)
            dd = dirs{d};
            legend_labels = cell(1,numel(algos));
            fig = figure;
            hold on;
            for k = 1:numel(algos)
                name = algos{k};
                bsd_e      = mean(e.bsd.(dd).(opt_curr_n).(name),  3);  % average over ears
                bsd_e_mu   = mean(bsd_e,2);
                bsd_e_std  = std(bsd_e, 0, 2);

                % Shaded error region
                    fill([fc; flipud(fc)], [bsd_e_mu - bsd_e_std; flipud(bsd_e_mu + bsd_e_std)], ...
                         colors(k,:), 'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
    
                % Plot solid line for Within FoV
                plot(fc, bsd_e_mu, '-o', 'LineWidth', 1.5, 'Color', colors(k,:));
                legend_labels{k} = name;
        
                
               
            end
            set(gca, 'XScale', 'log');
            legend(legend_labels, 'Location', 'northwest');
            xlim([20, 20e3]);
            ylim([0,10]);
            xlabel('Frequency (Hz)');
            ylabel('BSD (dB)');
            title(sprintf('BSD - %s of FoV', dirName{d}));
            grid on;
        
            % Save
            fig_name = fullfile(save_folder, sprintf('%s_%s_%s.fig', title_name, opt_curr_n , dd));
            savefig(fig_name);
            %close(fig);
        end
    end
  
end

function plot_errors_freq(BSMobj, e, title_name, save_folder)
%PLOT_ERRORS_FREQ Plots frequency-dependent error curves (NMSE & Normalized Magnitude) for four binaural algorithms
%   Creates 4 figures: NMSE Within FoV, NMSE Outside FoV, NMag Within FoV, NMag Outside FoV
%
%   Inputs:
%     - BSMobj: object with frequency vector in .freqs_sig
%     - e: struct with fields e.{nmse,nmag}.{BSM,WBSM,TBSM,CBSM}.{in,out}.mu/std (size [nFFT x 2])
%     - title_name: base title (string)
%     - save_folder: folder to save .fig files (string)

    algos    = fieldnames(e.nmse.in.ls);
    metrics  = {'nmse','mag'};
    dirs     = {'in','out'};
    optimization = {'ls','mls'};
    dirName  = {'Within','Outside'};

    % Frequency axis
    f = BSMobj.freqs_sig.';

    for opt_curr = 1:numel(optimization)
        opt_curr_n = optimization{opt_curr};
        for m = 1:numel(metrics)
            metric = metrics{m};
            for d = 1:numel(dirs)
                dd = dirs{d};
                fig = figure;
                hold on;
                legend_labels = cell(1,numel(algos));
                colors = lines(numel(algos));
    
                for k = 1:numel(algos)
                    name = algos{k};
                    % Extract per-ear data [nFFT x 2]
                    mu_ears  = e.(metric).(dd).(opt_curr_n).(name).mu;
                    std_ears = e.(metric).(dd).(opt_curr_n).(name).std;
                    % Average across ears
                    mu_lin  = mean(mu_ears, 2);
                    std_lin = mean(std_ears,2);
    
                    % Convert to dB
                    mu_db    = 10*log10(mu_lin);
                    sigma_db = 10*log10(mu_lin + std_lin) - mu_db;
    
                    % Shaded error region
                    fill([f; flipud(f)], [mu_db - sigma_db; flipud(mu_db + sigma_db)], ...
                         colors(k,:), 'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
    
                    % Mean curve
                    plot(f, mu_db, 'LineWidth',1.5,'Color',colors(k,:));
                    legend_labels{k} = name;
                end
    
                set(gca,'XScale','log');
                legend(legend_labels,'Location','northwest');
                xlim([20,20e3]);
                ylim([-40,12]);
                xlabel('Frequency (Hz)');
                ylabel(sprintf('%s (dB)', upper(metric)));
                title(sprintf('%s - %s of FoV', upper(metric), dirName{d}));
                grid on;
    
                % Save
                fig_name = fullfile(save_folder, sprintf('%s_%s_%s_%s.fig', title_name, opt_curr_n ,upper(metric), dd));
                savefig(fig, fig_name);
                %close(fig);
            end
        end
    end
end

function plot_ITD_curves(e, omega_az, title_name, save_folder)
    algos  = fieldnames(e);
    colors = lines(numel(algos));

    figure;
    hold on;
    legend_labels = {};

    for k = 1:numel(algos)
        name = algos{k};
        mu_in  = e.(name);  % average over ears

        plot(omega_az, mu_in, 'LineWidth', 1.5, 'Color', colors(k,:));
        legend_labels{end+1} = sprintf('%s', name);
    end

    legend(legend_labels, 'Location', 'northwest');
    xlim([omega_az(1),omega_az(end)]);
    ylim([-850,850]);
    xlabel("Incedent angle [deg]")
    ylabel("ITD [micro sec]")
    title(title_name);
    grid on;

    % Save
    fig_name = save_folder + title_name+".fig";
    savefig(fig_name);
end

function plot_ILD_curves(e, omega_az, title_name, save_folder)
    algos = fieldnames(e);
    colors = lines(numel(algos));

    figure;
    hold on;
    legend_labels = {};

    for k = 1:numel(algos)
        name = algos{k};
        mu_in  = squeeze(mean(e.(name),1));  % average over ears

        % Plot solid line for Within FoV
        plot(omega_az, mu_in, 'LineWidth', 1.5, 'Color', colors(k,:));
        legend_labels{end+1} = sprintf('%s', name);
    end

    legend(legend_labels, 'Location', 'northwest');
    xlim([omega_az(1),omega_az(end)]);
    ylim([-20,20]);
    xlabel("Incedent angle [deg]")
    ylabel("ILD [dB]")
    title(title_name);
    grid on;

    % Save
    fig_name = save_folder + title_name+".fig";
    savefig(fig_name);
end

function plot_ILD_error_curves(e, omega_az,fc, title_name, save_folder)
    algos = fieldnames(e);
    %algos(5) =[]; % remove JBSM
    %algos(3) =[]; % remove WBSM
    algos(1) =[]; % remove REF

    colors = lines(numel(algos));

    ref = e.REF;



    for k = 1:numel(algos)
        figure;
        name = algos{k};
        error = abs(ref - e.(name));
        error_smooth = imgaussfilt(error, 1.0); % requires Image Processing Toolbox
        imagesc(omega_az, fc, error_smooth);
        axis xy;                % Flip Y so that y_vec increases upwards
        xlabel('Incedent angle [deg]');
        ylabel('Frequncy [kHz]');
        c = colorbar;
        c.Label.String = 'Error [dB]';
        clim([0 5]);
        ylim([fc(1), fc(end)])
        xlim([omega_az(1), omega_az(end)])
        yticks([2000 5000 8000 10000 15000]);
        yticklabels(string(yticks));
        xticks([0 30 60 90 120 150 180]);
        xticklabels(string(xticks));
        title_name_tmp = title_name + "ILD error " + string(name);
        title(title_name_tmp);
        % Save
        fig_name = save_folder + title_name_tmp+".fig";
        savefig(fig_name);

    end



    
end

function res_out = package_results(objs, data, coeffs, errors,hrtf_hat)
%PACKAGE_RESULTS Collect outputs into a struct

    res_out.BSMobj = objs.BSM;
    res_out.c_BSM  = coeffs.BSM;
    res_out.c_WBSM = coeffs.WBSM;
    res_out.c_TBSM = coeffs.TBSM;
    res_out.c_JBSM = coeffs.JBSM;
    
    res_out.errors = errors;
    
    res_out.V_array_f_BSM = permute(data.V_k_bsm_normal,[1,3,2]);
    res_out.V_array_f_WBSM = permute(data.V_k_bsm,[1,3,2]);
    res_out.omega_array_BSM = data.omega_bsm_normal;
    res_out.omega_array_WBSM = data.omega_bsm;

    res_out.H_t = hrtf_hat.H_t;
    res_out.H_f = hrtf_hat.H_f;
end


function errors = evaluate_all_dirs_errors(errors, hrtf)
% EVALUATE_ALL_DIRS_ERRORS 
% Computes the global NMSE and Mag Error across ALL directions by 
% concatenating the Within-FoV and Outside-FoV predictions.

    disp("Compute Global (All Directions) NMSE and Magnitude Error")
    
    % Reconstruct the global reference HRTF
    ref_all = cat(1, hrtf.lebedev.in.ref, hrtf.lebedev.out.ref);
    
    % Define the algorithms and optimization types present in your struct
    algos = {'bsm', 'wbsm', 'tbsm', 'jbsm'};
    opts = {'ls', 'mls'};
    
    for a = 1:length(algos)
        algo = algos{a};
        algo_upper = upper(algo);
        
        for o = 1:length(opts)
            opt = opts{o};
            
            % Safely try to compute for each algorithm if it exists
            try
                % Concatenate the 'in' and 'out' predictions
                est_in  = hrtf.lebedev.in.(algo).(opt);
                est_out = hrtf.lebedev.out.(algo).(opt);
                est_all = cat(1, est_in, est_out);
                
                % Compute and store Global NMSE
                [nmse_mu, nmse_std] = calc_nmse(est_all, ref_all);
                errors.nmse.all.(opt).(algo_upper).mu  = nmse_mu;
                errors.nmse.all.(opt).(algo_upper).std = nmse_std;
                
                % Compute and store Global Magnitude Error
                [mag_mu, mag_std] = calc_mag(est_all, ref_all);
                errors.mag.all.(opt).(algo_upper).mu  = mag_mu;
                errors.mag.all.(opt).(algo_upper).std = mag_std;
            catch
                % Silently skip if the algorithm variant wasn't generated
            end
        end
    end
end


function plot_all_dirs_errors(BSMobj, e, title_name, save_folder, algos_to_plot)
% PLOT_ALL_DIRS_ERRORS Plots the global (All Directions) NMSE and 
% Magnitude errors. Creates a single figure 2x2 layout: 
% Rows = NMSE/Mag, Cols = LS/MLS.

    if ~isfield(e.nmse, 'all')
        warning('No "all" directions data found. Please run evaluate_all_dirs_errors first.');
        return;
    end

    % Default to plotting all available algorithms if none are specified
    if nargin < 5 || isempty(algos_to_plot)
        algos_to_plot = fieldnames(e.nmse.all.ls); 
    end

    f_sig = BSMobj.freqs_sig.'; % Frequency axis
    
    optimization = {'ls', 'mls'};
    opt_titles = {'LS', 'MagLS'};
    metrics = {'nmse', 'mag'};
    y_labels = {'NMSE (dB)', 'MagError (dB)'};

    fig = figure('Name', 'All Dirs Errors - LS vs MagLS'); 
    t = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % Generate a consistent color palette for the selected algorithms
    colors = lines(numel(algos_to_plot));
    
    for m = 1:numel(metrics)
        metric = metrics{m};
        
        for opt_idx = 1:numel(optimization)
            opt_type = optimization{opt_idx};
            
            ax = nexttile;
            hold on;
            
            % Check if this optimization type exists in the struct
            if ~isfield(e.(metric).all, opt_type)
                title(sprintf('%s - %s (No Data)', upper(metric), opt_titles{opt_idx}));
                continue;
            end
            
            for k = 1:numel(algos_to_plot)
                name = algos_to_plot{k};
                
                % Skip if this specific algorithm variant wasn't generated
                if ~isfield(e.(metric).all.(opt_type), name)
                    continue;
                end
                
                % Extract mean and std (averaging across the 2 ears)
                mu_lin  = mean(e.(metric).all.(opt_type).(name).mu, 2);
                std_lin = mean(e.(metric).all.(opt_type).(name).std, 2);
                
                % Convert to dB
                mu_plot = 10*log10(mu_lin);
                sigma   = 10*log10(mu_lin + std_lin) - mu_plot;
                
                % Shaded Error Region
                fill([f_sig; flipud(f_sig)], [mu_plot - sigma; flipud(mu_plot + sigma)], ...
                     colors(k,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

                % Mean Line
                plot(f_sig, mu_plot, 'LineWidth', 1.5, 'Color', colors(k,:), 'DisplayName', name);
            end
            
            % Formatting
            set(ax, 'XScale', 'log');
            grid on;
            xlim([20, 20e3]);
            ylim([-40, 10]); 
            
            % Dynamic Titles and Labels
            if m == 1
                title(sprintf('NMSE - All Dirs (%s)', opt_titles{opt_idx}));
            else
                title(sprintf('Mag Error - All Dirs (%s)', opt_titles{opt_idx}));
                xlabel('Frequency (Hz)');
            end
            
            if opt_idx == 1
                ylabel(y_labels{m});
            end
        end
    end
    
    % Global Legend at the bottom
    lgd = legend(algos_to_plot, 'Orientation', 'horizontal', 'Location', 'southoutside');
    lgd.Layout.Tile = 'south';

    % Save the figure
    fig_name = fullfile(save_folder, sprintf('%s_all_dirs_combined.fig', title_name));
    savefig(fig, fig_name);
end