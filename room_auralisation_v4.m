function res_out = room_auralisation_v4(roomsetup,arraysetup,simType)

shutup        = roomsetup.shutup;
R             = roomsetup.R;
S             = roomsetup.S;
freq_bands    = roomsetup.freq_bands;
roomDim       = roomsetup.roomDim;
arrayPos      = roomsetup.arrayPos;
sourcePos     = roomsetup.sourcePos;
arrayType     = arraysetup.arrayType;
sig_gender    = roomsetup.sig_gender;
direct_flag   = roomsetup.direct_flag;
diffuese_flag = roomsetup.diffuese_flag;
specular_flag = roomsetup.specular_flag;
fs_out        = roomsetup.fs;
roomsim_fs    = roomsetup.roomsim_fs;
N_PW          = roomsetup.N_PW;
c = roomsetup.c;
M = arraysetup.M;
if isfield(roomsetup, 'src_gain')
    src_gain = roomsetup.src_gain;
else
    src_gain = 1.0; 
end

switch simType
    case 'MCRoomSim'
    Room = SetupRoom('Dim', roomDim, ...
                     'Freq', freq_bands, ...
                     'Absorption', (1- R.^2) * ones(6, 6), ...
                     'Scattering', S * ones(6, 6)); % Increased scattering slightly to help diffusion
    
    Sources   = AddSource('Location', sourcePos, 'Type', sig_gender);
    Receivers = AddReceiver('Location', arrayPos, 'Type', 'sphharm', 'MaxOrder', N_PW, 'ComplexSH', false);
    Options   = MCRoomSimOptions('Fs', roomsim_fs,'Verbose',~shutup,'SimDirect',direct_flag,'SimSpec',specular_flag,'SimDiff',diffuese_flag,'SoundSpeed',c);  % consider adding Duration
    RIR_SH_t  = RunMCRoomSim(Sources, Receivers, Room, Options);
    
    % Run only with the direct reflecation
    Options         = MCRoomSimOptions('Fs', roomsim_fs,'Verbose',~shutup,'SimDirect',true,'SimSpec',false,'SimDiff',false,'SoundSpeed',c);  % consider adding Duration
    RIR_SH_t_direct = RunMCRoomSim(Sources, Receivers, Room, Options);
    
    
    % Adjust SH orders according to our conventions:
    P = MCRoomPerm(N_PW)';
    C = SHc2r(N_PW)';
    RIR_SH_t = (C * P * RIR_SH_t.').';
    RIR_SH_t_direct = (C * P * RIR_SH_t_direct.').';

    DRR = RoomParams.DRR(RIR_SH_t(:,1), RIR_SH_t_direct(:,1), "outputType", "dB");
    T60 = RoomParams.T60(RIR_SH_t(:,1), roomsim_fs);
    Rc  = RoomParams.critical_distance_diffuse(roomDim, 1-R);
    if ~shutup
        PlotSimSetup(Sources,Receivers,Room)
    end



    clear P C
    case 'Tom'
        refCoef = R;
        N_INFTY = N_PW;
        Anechoic = ~(diffuese_flag || specular_flag);
        [RIR_SH_t, tmp1,tmp2] = image_method.calc_rir(roomsim_fs, roomDim, sourcePos, arrayPos, refCoef, {"max_reflection_order",50,"angle_dependence",true,"zero_first_delay",false}, {"array_type", "anm", "N", N_INFTY},Anechoic);

        DRR = tmp2.DRR;
        T60 = tmp2.T60;
        Rc  = tmp2.rd;
end

% --- NEW: Apply Source Gain Here ---
% Since the system is linear, this is mathematically perfect for Omni sources.
RIR_SH_t = RIR_SH_t * src_gain;



if fs_out ~= roomsim_fs
    gComDiv = gcd(roomsim_fs, fs_out);
    p = double(fs_out / gComDiv);
    q = double(roomsim_fs / gComDiv);
     % Resample the SH RIRs (Optional, but good for completeness)
     RIR_SH_t  = resample(RIR_SH_t, p, q);
end








if diffuese_flag || specular_flag
    t_vec = (1:size(RIR_SH_t,1))*(1/fs_out);
    [~,min_idx] = min(abs(t_vec - T60));
    RIR_SH_t = RIR_SH_t(1:min_idx,:);
end


% if ~shutup
%     figure
%     plot((0:size(RIR_SH_t,1)-1)/fs_out, real(RIR_SH_t(:,1))); 
%     xlabel('Time [sec]'); % plot the RIR of a00
% 
%     %PlotSimSetup(Sources, Receivers, Room)
% end

% Spherical coordinates of direct sound 
direct_sound_rel_cart = sourcePos - arrayPos;
[th0, ph0, r0]=c2s(direct_sound_rel_cart(1), direct_sound_rel_cart(2), direct_sound_rel_cart(3));
ph0 = mod(ph0, 2*pi);
direct_sound_rel_sph = [r0, th0, ph0];

if ~shutup
    disp(['Source position: (r,th,ph) = (' num2str(direct_sound_rel_sph(1),'%.2f') ','...
        num2str(direct_sound_rel_sph(2)*180/pi,'%.2f') ','...
        num2str(direct_sound_rel_sph(3)*180/pi,'%.2f') ')']);
        
    fprintf('MCRoomSim T60 = %.2f sec\n', T60);
    fprintf('MCRoomSim Critical Distance = %.2f m\n', Rc);
    fprintf('MCRoomSim DRR = %.2f ', DRR);
end

if ~shutup
    fprintf('Finished MCRoomSim RIR generation\n');
    fprintf('Start calculating array measurements (frequency domain)\n');
end
f_cut_magLS = 1.5e3;
at = arraysetup.arrayType;
az = 0;
nfft = 512;
filt_len = nfft / fs_out;
angles = 0;
BSMobj           = setup_BSM_object(M,at,az,fs_out,filt_len,4,N_PW,f_cut_magLS,angles);
omega_bsm        = [BSMobj.th_BSMgrid_vec,BSMobj.ph_BSMgrid_vec,ones(size(BSMobj.ph_BSMgrid_vec))];


if arrayType == 5
    load([BSMobj.atf_folder,BSMobj.atf_name]);
    th_atf = ATF.th;
    ph_atf = ATF.ph;
    Omega_bsm_2 = [th_atf,mod(ph_atf,2*pi)];
end


[Vf, Vt] = get_me_some_V(BSMobj,omega_bsm);
Vf       = permute(Vf,[1,3,2]);

NFFT_IR = 2^nextpow2(size(RIR_SH_t,1));
NFFT_V = (size(Vf,2)-1)*2;

NFFT = max(NFFT_IR,NFFT_V);
RIR_SH_f = fft(RIR_SH_t, NFFT, 1);
% remove negative frequencies
RIR_SH_f(NFFT/2+2:end,:)=[];
nFreq = size(RIR_SH_f, 1);          % 

f_RIR = linspace(0, fs_out/2, nFreq).';       % Frequencies of RIR_SH_f
f_Vf = BSMobj.freqs_sig.';                        % Frequencies at which Vf is defined (e.g. 257 values from measured HRTFs)
f_RIR(1) = f_Vf(1);                               % Frequencies of RIR_SH_f


I = BSMobj.M;

% --- FIX THE LEFT/RIGHT MIRROR ---
% Invert the Azimuth to align the SH basis with the Room Simulator's handedness
ph_mirrored = 2*pi - omega_bsm(:,2); 
Y_a = sh2(N_PW, omega_bsm(:,1), ph_mirrored);
% ---------------------------------

Yp_a = pinv(Y_a);  % size (Q x (N+1)^2)

nSH = (N_PW+1)^2;
Q = size(Vf,3);


% % Vf: [I x length(f_Vf) x Q]
% Vf_interp = zeros(I, nFreq, Q);
% for i = 1:I
%     for q = 1:Q
%         %Vf_interp(i,:,q) = interp1(f_Vf, squeeze(Vf(i,:,q)), f_RIR, 'linear', 'extrap');
%         Vf_interp(i,:,q) = interp1(f_Vf, squeeze(Vf(i,:,q)), f_RIR, 'pchip');
%         %interp1(f_Vf, Vf, f_RIR, 'pchip');
%     end
% end
% Vf = Vf_interp ;

% --- EXACT FFT INTERPOLATION (Replaces PCHIP) ---
% 1. Reconstruct full spectrum of the ATF
NFFT_V = (size(Vf,2)-1)*2;
Vf_full = cat(2, Vf, conj(Vf(:, end-1:-1:2, :)));

% 2. Convert back to Time Domain [I x NFFT_V x Q]
vt_time = real(ifft(Vf_full, NFFT_V, 2)); 
%vt_time = circshift(vt_time,30,2);
vt_time = circshift(vt_time,NFFT_V/2,2);
%vt_time = fftshift(vt_time,2);

% 3. Pad with zeros to match the huge RIR NFFT
vt_padded = zeros(I, NFFT, Q);
vt_padded(:, 1:NFFT_V, :) = vt_time;
vt_padded = circshift(vt_padded,-1*NFFT_V/2,2);

% 4. Convert back to Frequency Domain at high resolution
Vf_highres = fft(vt_padded, NFFT, 2);

% 5. Extract positive frequencies to match RIR_SH_f
Vf_interp = Vf_highres(:, 1:nFreq, :);
Vf = Vf_interp;
% ------------------------------------------------






p_f = zeros(nFreq, I);  % Frequency domain pressure measurements (Freq x Mics)

for i = 1:I
    % Compute Hnm_f (Freq x SH components)
    Hnm_f = squeeze(Vf(i,:,:)) * Yp_a;  % size: Freq x nSH
    for f = 1:nFreq
        p_f(f,i) = Hnm_f(f,:) * RIR_SH_f(f,:).';
    end
end



% Check the array time doamin responce
p_array_f = p_f;
p_array_f(end+1:NFFT, :) = 0;
p_array_t = ifft(p_array_f, "symmetric");
p_array_f = p_f; %use only positive bins


% if fs_out ~= roomsim_fs
%     if ~shutup
%         fprintf('Resampling final output from %d Hz to %d Hz...\n', sim_fs, target_fs);
%     end
% 
%     [p, q] = rat(fs_out / roomsim_fs);
% 
%     % Resample the Microphone Signals
%     p_array_t = resample(p_array_t, p, q);
% 
%     % Resample the SH RIRs (Optional, but good for completeness)
%     RIR_SH_t  = resample(RIR_SH_t, p, q);
% 
%     % Re-calculate Frequency Domain (to match the new time signal)
%     % This ensures p_array_f corresponds to the 16k signal, not the 48k one
%     NFFT_new  = 2^nextpow2(size(p_array_t,1));
%     p_array_f = fft(p_array_t, NFFT_new);
%     p_array_f(NFFT_new/2+2:end, :) = []; % Keep positive freq only
% end


res_out.p_array_t = p_array_t;
res_out.p_array_f = p_array_f;
res_out.RIR_SH_t = RIR_SH_t;
res_out.T60 = T60;
res_out.T60_target = roomsetup.T60_target;
res_out.Rc = Rc;
res_out.DRR = DRR;
res_out.src_gain  = src_gain;
res_out.fs = fs_out;

if ~shutup
    fprintf('Finished calculating array measurements (frequency domain)\n');
end