function BSMobj= setup_BSM_object(M,arrayType,head_rot_az,fs,filt_len,source_distribution,Q,f_cut_magLS,sector_ang)

%
% filt_len,source_distribution,Q,N_PW
% parameters/flags - array
% filt_len = 0.032;                                      % filters (BSM/HRTF) length [sec]
% filt_len = filt_len/2;

rigidArray = 1;                                        % 0 - open array, 1 - rigid array
%r_array = 0.1;                                         % array radius

%normSV = true;                                         % true - normalize steering vectors
normSV = true;                                         % true - normalize steering vectors

% parameters/flags - general
c = 343;                                               % speed of sound [m/s]
%fs = 44100;                                           % choose samplong frequency in Hz
% N_PW = 14;                                             % SH order of plane-wave synthesis %HUTUBS=16, META=13, SONICOM=7. Check condition number of Y to set this parameter, it dependes on the space sample schem of the HRTF
% if arrayType == 3
%     % This is for EasyCom Glass
%     N_PW = 14;
% end

% parameters/flags - BSM design
inv_opt = 1;                                           % opt=1 -> ( (1 / lambda) * (A * A') + eye )  |||| opt2=1 -> ((A * A') + lambda * eye);
%source_distribution = 2;                               % 0 - nearly uniform (t-design), 1 - spiral nearly uniform
%Q = N_PW;                                               % Assumed number of sources

% parameters/flags - MLS solver
tol_magLS = 1e-20;    
max_iter_magLS = 1E5;

%noise
SNR = 20;                                              % assumed sensors SNR [dB]
sig_n = 1;
sig_s = 10^(SNR/10) * sig_n;
SNR_lin = sig_s / sig_n;    

%signal
filt_samp    = ceil(filt_len * fs);
freqs_sig    = ( 0 : (filt_samp / 2) ) * fs / filt_samp;
freqs_sig(1) = 1/4 * freqs_sig(2); %to not divide by zero



% Text variables for plots 
if ~rigidArray
    sphereType = "open";
else
    sphereType = "rigid";
end
switch arrayType 
    case 0
        arrayTypeTxt = sphereType + "Spherical";
        %N_PW = 18;
    case 1
        arrayTypeTxt = sphereType + "SemiCirc";
        %N_PW = 18;
    case 2
        arrayTypeTxt = sphereType + "FullCirc";
        %N_PW = 18;
    case 3
        arrayTypeTxt = "Glasses";
        atf_folder        = '/Users/orberebi/Documents/GitHub/PBSM/ATFs/';
        atf_name          = 'Device_ATFs_v2.h5';
        BSMobj.atf_folder = atf_folder;
        BSMobj.atf_name   = atf_name;
        %N_PW = 18;
        M = 4;
    case 4
        arrayTypeTxt = sphereType + "Sim_Glass";
        atf_folder        = '/Users/orberebi/Documents/GitHub/PBSM/ATFs/';
        atf_name          = 'Device_ATFs_v2.h5';
        BSMobj.atf_folder = atf_folder;
        BSMobj.atf_name   = atf_name;
        %N_PW = 18;
        M = 4;
    case 5
        arrayTypeTxt = sphereType + "PVT";
        atf_folder        = '/Users/orberebi/Documents/GitHub/PBSM/ATFs/';
        atf_name          = 'SN_PVT_HATS_atf_12-8-23.mat';
        BSMobj.atf_folder = atf_folder;
        BSMobj.atf_name   = atf_name;
        %N_PW = 18;
        M = 5;
    case 6 % Aria
        arrayTypeTxt = "Aria";
        atf_folder        = '/Users/orberebi/Documents/GitHub/IWAENC_2026_WBSM/ATF/';
        if M == 5
            atf_name          = 'aria_interpolated_M5.mat';
        else
            atf_name          = 'aria_interpolated.mat';
            M = 7;
        end
        BSMobj.atf_folder = atf_folder;
        BSMobj.atf_name   = atf_name;
        
        
        %M = 5;

end

if fs > 16e3
    N_PW = 35;
else
    N_PW = 18;
end

%disp("Num of microphones = "+ num2str(M))               % number of microphones
%disp("Mic Array type = " + arrayTypeTxt)                % 0 - spherical array, 1 - semi-circular array, 2 - full-circular array
%disp("HRTF rotation = " + num2str(rad2deg(head_rot_az)))

[th_BSMgrid_vec, ph_BSMgrid_vec] = BSM_toolbox.BSMgrid(source_distribution, Q,sector_ang);
r_BSMgrid_vec = ones(size(th_BSMgrid_vec));
% Get a new BSMgrid 



n_mic = M;
Q = length(th_BSMgrid_vec);
%[th_array, ph_array, ~] = BSM_toolbox.GetArrayPositions(arrayType, n_mic, 0);       
[th_array, ~, ph_array,r_array] = BSM_toolbox.GetArrayPositions(arrayType, n_mic, head_rot_az);


BSMobj.freqs_sig            = freqs_sig;
BSMobj.head_rot_az          = head_rot_az;
BSMobj.N_PW                 = N_PW;    
BSMobj.c                    = c;
BSMobj.r_array              = r_array;
BSMobj.rigidArray           = rigidArray;
BSMobj.th_BSMgrid_vec       = th_BSMgrid_vec;
BSMobj.ph_BSMgrid_vec       = ph_BSMgrid_vec;
BSMobj.r_BSMgrid_vec        = r_BSMgrid_vec;
BSMobj.Tik                  = [];
BSMobj.Jan                  = [];
BSMobj.beta                 = [];
BSMobj.sector_ang           = sector_ang;
BSMobj.f_cut_magLS          = f_cut_magLS;
BSMobj.tol_magLS            = tol_magLS;
BSMobj.max_iter_magLS       = max_iter_magLS;
BSMobj.normSV               = normSV;
BSMobj.SNR_lin              = SNR_lin;
BSMobj.inv_opt              = inv_opt;
BSMobj.M                    = M;
BSMobj.Q                    = Q;
BSMobj.source_distribution  = source_distribution;
BSMobj.fs                   = fs;
BSMobj.filt_samp            = filt_samp;
BSMobj.sphereType           = sphereType;
BSMobj.arrayTypeTxt         = arrayTypeTxt;
BSMobj.arrayType            = arrayType;
BSMobj.n_mic                = n_mic;
BSMobj.th_array             = th_array;
BSMobj.ph_array             = ph_array;
BSMobj.desired_fs           = fs;
BSMobj.magLS_cvx            = false;

end
