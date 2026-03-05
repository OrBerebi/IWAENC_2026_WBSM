clear
clc 
close all
restoredefaultpath;

base_path =  string(pwd);
startup_script(base_path)


% Define the *relative* path to the results folder within base_path
figures_save_folder = "/results/26_02_26_yonatan_stuff/";

%run_beta_ablation(base_path, figures_save_folder);

% WBSM parameter
beta = 0.075;

% array config
arraysetup.M            = 6;
arraysetup.arrayType    = 1;
arraysetup.sphereType   = "rigid";
arraysetup.fs           = 48e3; 
[arraysetup.th_array, ~, arraysetup.ph_array, arraysetup.r_array] = ...
    BSM_toolbox.GetArrayPositions(arraysetup.arrayType, arraysetup.M, 0);

fprintf('Pre-loading data structures...\n');
dummy_params    = init_params_local(1, base_path, figures_save_folder);
dummy_params.fs = arraysetup.fs;
data_static     = load_and_prepare_data_local(dummy_params); 

plot_anechoic_analysis = true;
res_out  = anechoic_analysis_FOV_BSM_v5(beta,arraysetup, base_path, figures_save_folder, data_static, plot_anechoic_analysis);





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
    % Azimuth plane grid
    az_res = 1;
    %phi_az = deg2rad(linspace(0,180,180/az_res).');
    phi_az = deg2rad(linspace(-180,180,360/az_res).');

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

