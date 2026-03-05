function V_k = CalculateSteeringVectors(BSMobj, N_SV, th_DOA, ph_DOA)
%% CalculateSteeringVectors.m
% Calculate steering vectors of supported arrays (spherical / semi-circular / fully-circular)
% O.B added support for easy-com Glasses ATF
% Equation for calculation is: Ymic * B * Y'
%% Inputs:
% BSMobj            : (MATLAB object) containing parameters
% N_SV              : (scalar) maximal SH order of calculated steering vectors
% th_DOA            : (DOAs x 1) elevation of arriving/incident plane waves [rad]
% ph_DOA            : (DOAs x 1) azimuth of arriving/incident plane waves   [rad]
%% Outputs:
% V_k               : (n_mic x Q x freq) steering vectors (ATF)
% V_k               : (freq x DOAs x M) steering vectors in freq. domain

    %init
    r_array         = BSMobj.r_array;    
    freqs_sig       = BSMobj.freqs_sig;    
    c               = BSMobj.c;
    n_mic           = BSMobj.n_mic;                      
    normSV          = BSMobj.normSV;
    sphereType      = BSMobj.sphereType;    
    
    if strcmp(BSMobj.arrayTypeTxt, "Glasses")
        atf_folder = BSMobj.atf_folder;
        atf_name   = BSMobj.atf_name;
        th_atf = h5read(fullfile(atf_folder,atf_name),'/Theta');
        ph_atf = h5read(fullfile(atf_folder,atf_name),'/Phi');
        Vt = h5read(fullfile(atf_folder,atf_name),'/IR');
        Vt = permute(Vt,[1,3,2]);
        atf_fs = h5read(fullfile(atf_folder,atf_name),'/SamplingFreq_Hz');
        if atf_fs ~= BSMobj.fs
            Vt = resample(Vt, BSMobj.fs, atf_fs, Dimension=2);
        end
        %Mics 5 and 6 are inside the ears - thus, remove them from simulation
        %Vt = Vt(1:4, :, :);
        if BSMobj.n_mic == 4
            Vt = Vt(1:4, :, :);
        end
        M = size(Vt, 1);
        BSMobj.n_mic = M;
    
        % ================= ATFs preprocessing ============
        % Interpolate ATF to BSMGrid
        
        % Convert ATF to SH domain
        N_ATF = 16;
        % Y_atf_pinv = pinv(shmat(N_ATF, [th_atf, ph_atf]));
        Y_atf_before_pinv = shmat(N_ATF, [th_atf, ph_atf]);
        
        [Vt_m, Vt_s, Vt_d] = size(Vt);
        Vt_nm = reshape(Vt, [Vt_m * Vt_s, Vt_d]);
        % Vt_nm = (Y_atf_pinv * Vt_nm.').';
        Vt_nm = (Y_atf_before_pinv \ (Vt_nm.')).';
        
        % Convert back to space domain in DOA grid
        Y_atf = shmat(N_ATF, [th_DOA, ph_DOA]);
        Vt_BSM = (Y_atf * Vt_nm.').';
        Vt_BSM = reshape(Vt_BSM, [Vt_m, Vt_s, length(th_DOA)]);
        % Vt_BSM is [f_len x directions x n_mic]
    
        % Convert ATF to frequency domain
        NFFT = BSMobj.filt_samp;
        Vf_BSM = fft(Vt_BSM, NFFT, 2);
        
        % remove negative frequencies
        Vf_BSM(:, NFFT/2+1 + 1:end, :)=[]; % error of lior? +1 forgot
        %freqs_atf = (0 : (size(Vf_BSM, 2) - 1))'*(BSMobj.fs / NFFT);
        
        % steering vector should be [n_mic x directions x freq] so reshape
        V_k = permute(Vf_BSM, [2, 3, 1]);
    else
        
        th_array        = BSMobj.th_array;
        ph_array        = BSMobj.ph_array;  
        if size(th_DOA, 2) > size(th_DOA, 1)
            th_DOA = th_DOA .';
        end    
        if size(ph_DOA, 2) > size(ph_DOA, 1)
            ph_DOA = ph_DOA .';
        end
    
        f_len = length(freqs_sig);
        N_len = (N_SV + 1)^2;
        DOAs_len = length(th_DOA);
        Y = conj(shmat(N_SV, [th_DOA, ph_DOA], true, true));            % [SH x PWs]
        Ymic   = shmat(N_SV, [th_array.', ph_array.'], true, false);    % [mic x SH]  
        
        kr = 2 * pi * r_array * freqs_sig / c;
        b = bn(N_SV, kr, "sphereType", sphereType);            
        bb = repmat(b, 1, 1, DOAs_len);
        YY = repmat(Y, 1, 1, f_len); YY = permute(YY, [3 1 2]);
        bY = bb .* YY;
        bY1 = permute(bY, [2 1 3]); 
        bY2 = reshape(bY1, N_len, f_len * DOAs_len);
        V_k = Ymic * bY2;
        V_k = reshape(V_k, n_mic, f_len, DOAs_len);
        V_k = permute(V_k, [2 3 1]);
        % steering vector is [f_len x directions x n_mic]
        
    
    end
    %normalize steering vectors    
    if normSV
        V_k_norm = vecnorm(V_k, 2, 3);
        V_k = V_k ./ V_k_norm;        
    end
    
end