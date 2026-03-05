function [V_k, V_t] = CalculateSteeringVectors_v2(BSMobj, N_SV, th_DOA, ph_DOA)
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
    Vt_BSM = [];

    if strcmp(BSMobj.arrayTypeTxt, "Aria")
        nFFT = BSMobj.filt_samp;
        %aria_atf = sofaToEaro(convertStringsToChars([BSMobj.atf_folder,BSMobj.atf_name]));
        aria_atf = load(string([BSMobj.atf_folder,BSMobj.atf_name])).aria;
        %aria_lebedev = interpolate_aria(aria_atf);
        %aria_lebedev.data = aria_lebedev.data(:,1:2^8,:);
        %aria_lebedev.taps = 2^8;
        Vt = aria_atf.IR_interpolated;
        if aria_atf.fs ~= BSMobj.fs
            gComDiv = gcd(aria_atf.fs, BSMobj.fs);
            p = double(BSMobj.fs / gComDiv);
            q = double(aria_atf.fs / gComDiv);
            % Resample the SH RIRs (Optional, but good for completeness)
            newData = zeros(size(Vt,1), ceil(size(Vt,2) * (p/q)), size(Vt,3));
            for ii=1:size(Vt,1)
               newData(ii,:,:) = resample(squeeze(Vt(ii,:,:)),p,q);
            end
            Vt = newData;
        end


        aria_atf.IR = Vt;
        aria_atf.IR = cat(2, aria_atf.IR, zeros(size(aria_atf.IR,1),nFFT-size(aria_atf.IR,2),size(aria_atf.IR,3)));
        Vt = aria_atf.IR;
        th_atf = aria_atf.lebedev_elevation;
        ph_atf = aria_atf.lebedev_azimuth;
        
        Vt = permute(Vt,[3,2,1]); %Mics x Time samples x Directions





        % ================= ATFs preprocessing ============
        % Interpolate ATF to BSMGrid
        % Convert ATF to SH domain
        N_ATF = 32;
        %Y_atf_pinv = pinv(shmat(N_ATF, [th_atf, ph_atf]));
        Yq = pinv(sh2(N_ATF, th_atf, ph_atf));
        Y_atf_before_pinv = shmat(N_ATF, [th_atf, ph_atf]);
        
        [Vt_m, Vt_s, Vt_d] = size(Vt);
        Vt_nm = reshape(Vt, [Vt_m * Vt_s, Vt_d]);
        Vt_nm = (Y_atf_before_pinv \ (Vt_nm.')).';

        for nn=1:Vt_m       % iterate through mics
            tmp(nn,:,:) = (double(squeeze(Vt(nn,:,:)))*Yq);
        end
        Vt_nm_new = tmp; % mics x time x SH
        clear tmp
        % Convert back to space domain in DOA grid
        if size(th_DOA, 2) > size(th_DOA, 1)
            th_DOA = th_DOA .';
        end    
        if size(ph_DOA, 2) > size(ph_DOA, 1)
            ph_DOA = ph_DOA .';
        end
        
        Y_atf = shmat(N_ATF, [th_DOA, ph_DOA]);
        Vt_BSM = (Y_atf * Vt_nm.').';
        Vt_BSM = reshape(Vt_BSM, [Vt_m, Vt_s, length(th_DOA)]);
        Vt_BSM = real(Vt_BSM);


        Y = sh2(N_ATF, th_DOA, ph_DOA);
         % Iterate through sources
        for nn=1:Vt_m
            tmp(nn,:,:) = double(squeeze(Vt_nm_new(nn,:,:))*Y);
        end
        Vt_BSM_new = tmp;
        clear tmp



        %Switch between verstions
        %Vt_BSM = Vt_BSM_new;

        % Vt_BSM is [f_len x directions x n_mic]
    
        % Convert ATF to frequency domain
        NFFT = BSMobj.filt_samp;
        Vf_BSM = fft(Vt_BSM, NFFT, 2);
        
        % remove negative frequencies
        Vf_BSM(:, NFFT/2+1 + 1:end, :)=[]; % error of lior? +1 forgot
        %freqs_atf = (0 : (size(Vf_BSM, 2) - 1))'*(BSMobj.fs / NFFT);
        
        % steering vector should be [n_mic x directions x freq] so reshape
        V_k = permute(Vf_BSM, [2, 3, 1]);



    elseif strcmp(BSMobj.arrayTypeTxt, "Glasses")
        atf_folder = BSMobj.atf_folder;
        atf_name   = BSMobj.atf_name;
        th_atf = h5read(fullfile(atf_folder,atf_name),'/Theta');
        ph_atf = h5read(fullfile(atf_folder,atf_name),'/Phi');
        Vt = h5read(fullfile(atf_folder,atf_name),'/IR');
        Vt = permute(Vt,[1,3,2]); %Mics x Time samples x Directions
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
        N_ATF = N_SV;
        %Y_atf_pinv = pinv(shmat(N_ATF, [th_atf, ph_atf]));
        Yq = pinv(sh2(N_ATF, th_atf, ph_atf));
        Y_atf_before_pinv = shmat(N_ATF, [th_atf, ph_atf]);
        
        [Vt_m, Vt_s, Vt_d] = size(Vt);
        Vt_nm = reshape(Vt, [Vt_m * Vt_s, Vt_d]);
        %Vt_nm = (Y_atf_pinv * Vt_nm.').';
        Vt_nm = (Y_atf_before_pinv \ (Vt_nm.')).';

        for nn=1:Vt_m       % iterate through mics
            tmp(nn,:,:) = (double(squeeze(Vt(nn,:,:)))*Yq);
        end
        Vt_nm_new = tmp; % mics x time x SH
        clear tmp
        % Convert back to space domain in DOA grid
        if size(th_DOA, 2) > size(th_DOA, 1)
            th_DOA = th_DOA .';
        end    
        if size(ph_DOA, 2) > size(ph_DOA, 1)
            ph_DOA = ph_DOA .';
        end
        
        Y_atf = shmat(N_ATF, [th_DOA, ph_DOA]);
        Vt_BSM = (Y_atf * Vt_nm.').';
        Vt_BSM = reshape(Vt_BSM, [Vt_m, Vt_s, length(th_DOA)]);


        Y = sh2(N_ATF, th_DOA, ph_DOA);
         % Iterate through sources
        for nn=1:Vt_m
            tmp(nn,:,:) = double(squeeze(Vt_nm_new(nn,:,:))*Y);
        end
        Vt_BSM_new = tmp;
        clear tmp

        %Switch between verstions
        %Vt_BSM = Vt_BSM_new;

        % Vt_BSM is [f_len x directions x n_mic]
    
        % Convert ATF to frequency domain
        NFFT = BSMobj.filt_samp;
        Vf_BSM = fft(Vt_BSM, NFFT, 2);
        
        % remove negative frequencies
        Vf_BSM(:, NFFT/2+1 + 1:end, :)=[]; % error of lior? +1 forgot
        %freqs_atf = (0 : (size(Vf_BSM, 2) - 1))'*(BSMobj.fs / NFFT);
        
        % steering vector should be [n_mic x directions x freq] so reshape
        V_k = permute(Vf_BSM, [2, 3, 1]);
    elseif strcmp(BSMobj.arrayTypeTxt, "rigidPVT")
        atf_folder = BSMobj.atf_folder;
        atf_name   = BSMobj.atf_name;
        load([atf_folder,atf_name]);
        th_atf = ATF.th;
        ph_atf = ATF.ph;
        
        % Convert back to space domain in DOA grid
        if size(th_DOA, 2) > size(th_DOA, 1)
            th_DOA = th_DOA .';
        end    
        if size(ph_DOA, 2) > size(ph_DOA, 1)
            ph_DOA = ph_DOA .';
        end

        Vt = ATF.IR;
        Vt = permute(Vt,[1,3,2]); %Mics x Time samples x Directions
        atf_fs = ATF.params.samplerate;
        if atf_fs ~= BSMobj.fs
            Vt = resample(Vt, BSMobj.fs, atf_fs, Dimension=2);
        end
        M = size(Vt, 1);
        BSMobj.n_mic = M;
        % ================= ATFs preprocessing ============
        % Interpolate ATF to BSMGrid
        % Convert ATF to SH domain
        N_ATF = N_SV;
        
        %Y_atf_pinv = pinv(shmat(N_ATF, [th_atf, ph_atf]));
        Yq = pinv(sh2(N_ATF, th_atf, ph_atf));
        
        %Old method
        % [Vt_m, Vt_s, Vt_d] = size(Vt);
        % Y_atf_before_pinv = shmat(N_ATF, [th_atf, ph_atf]);
        % Vt_nm = reshape(Vt, [Vt_m * Vt_s, Vt_d]);
        % %Vt_nm = (Y_atf_pinv * Vt_nm.').';
        % Vt_nm = (Y_atf_before_pinv \ (Vt_nm.')).';


        for nn=1:M       % iterate through mics
            tmp(nn,:,:) = (double(squeeze(Vt(nn,:,:)))*Yq);
        end
        Vt_nm_new = tmp; % mics x time x SH
        clear tmp
        
        
        %Old method
        % Y_atf = shmat(N_ATF, [th_DOA, ph_DOA]);
        % Vt_BSM = real((Y_atf * Vt_nm.').');
        % Vt_BSM = reshape(Vt_BSM, [Vt_m, Vt_s, length(th_DOA)]);


        Y = sh2(N_ATF, th_DOA, ph_DOA);
         
        for nn=1:M  % iterate through mics
            tmp(nn,:,:) = double(squeeze(Vt_nm_new(nn,:,:))*Y);
        end
        Vt_BSM_new = tmp;
        clear tmp

        %Switch between verstions
        Vt_BSM = Vt_BSM_new;

        % Vt_BSM is [f_len x directions x n_mic]
    
        % Convert ATF to frequency domain
        NFFT = BSMobj.filt_samp;
        Vf_BSM = fft(Vt_BSM, NFFT, 2);
        
        % remove negative frequencies
        Vf_BSM(:, (NFFT/2)+1 + 1:end, :)=[]; % error of lior? +1 forgot
        %freqs_atf = (0 : (size(Vf_BSM, 2) - 1))'*(BSMobj.fs / NFFT);
        
        % steering vector should be [freq x directions x n_mic] so reshape
        V_k = permute(Vf_BSM, [2, 3, 1]);



    else
        %N_SV = 35;
        th_array        = BSMobj.th_array;
        ph_array        = BSMobj.ph_array;
        
        if size(th_DOA, 2) > size(th_DOA, 1)
            th_DOA = th_DOA .';
        end    
        if size(ph_DOA, 2) > size(ph_DOA, 1)
            ph_DOA = ph_DOA .';
        end
    
        f_len    = length(freqs_sig);
        N_len    = (N_SV + 1)^2;
        DOAs_len = length(th_DOA);
        Y        = conj(shmat(N_SV, [th_DOA, ph_DOA], true, true));       % [SH x PWs]
        Ymic     = shmat(N_SV, [th_array.', ph_array.'], true, false);    % [mic x SH]  
        
        kr  = 2 * pi * r_array(1) * freqs_sig / c;
        b   = bn(N_SV, kr, "sphereType", sphereType);            
        bb  = repmat(b, 1, 1, DOAs_len);
        YY  = repmat(Y, 1, 1, f_len); YY = permute(YY, [3 1 2]);
        bY  = bb .* YY;
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

    if nargout > 1
        Vt_BSM = permute(Vt_BSM,[2,3,1]);
        V_t = Vt_BSM;
    end
    
end