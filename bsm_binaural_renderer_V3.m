function res_out = bsm_binaural_renderer_V3(BSMobj,c_BSM_FOV,c_BSM,res_out_room,roomsetup)



fs        = BSMobj.fs;
f_Vf      = BSMobj.freqs_sig.';                % Frequencies at which Vf is defined (e.g. 257 values from measured HRTFs)
nFreq     = size(res_out_room.p_array_f, 1);
f_RIR     = linspace(0, fs/2, nFreq).';       % Frequencies of RIR_SH_f
f_RIR(1)  = f_Vf(1);                       % Frequencies of RIR_SH_f
NFFT      = size(res_out_room.p_array_t,1);
NFFT_V    = BSMobj.filt_samp;
%N_PW      = BSMobj.N_PW;
sig_path  = roomsetup.sig_path;
HRTF_path = roomsetup.HRTF_path;

I = BSMobj.M;
Q = size(BSMobj.th_BSMgrid_vec,1);

p_array_f = res_out_room.p_array_f;

RIR_SH_t = res_out_room.RIR_SH_t;


[s, signal_fs] = audioread(sig_path);
if (signal_fs ~= fs)
    s = resample(s, fs, signal_fs);
    signal_fs = fs;
end

%fprintf('Start calculating BRIRs\n');

c_ls_FOV = c_BSM_FOV.f_ls; c_ls_FOV = permute(c_ls_FOV,[2,1,3]);
c_mls_FOV = c_BSM_FOV.f_mls; c_mls_FOV = permute(c_mls_FOV,[2,1,3]);

c_ls = c_BSM.f_ls; c_ls = permute(c_ls,[2,1,3]);
c_mls = c_BSM.f_mls; c_mls = permute(c_mls,[2,1,3]);



% Interpulate the BSM coeff to the right NFFT bins
f_BSM = BSMobj.freqs_sig.';
c_interp_ls_f = zeros(nFreq, I, 2);
c_interp_mls_f = zeros(nFreq, I, 2);
c_interp_ls = zeros(nFreq, I, 2);
c_interp_mls = zeros(nFreq, I, 2);

for i = 1:I
    for q = 1:2
        c_interp_ls_f(:,i,q)    = interp1(f_BSM, squeeze(c_ls_FOV(:,i,q)), f_RIR, 'pchip');
        c_interp_mls_f(:,i,q)   = interp1(f_BSM, squeeze(c_mls_FOV(:,i,q)), f_RIR, 'pchip');
        c_interp_ls(:,i,q)      = interp1(f_BSM, squeeze(c_ls(:,i,q)), f_RIR, 'pchip');
        c_interp_mls(:,i,q)     = interp1(f_BSM, squeeze(c_mls(:,i,q)), f_RIR, 'pchip');

    end
end

c_ls_FOV  = c_interp_ls_f;
c_mls_FOV = c_interp_mls_f;
c_ls      = c_interp_ls;
c_mls     = c_interp_mls;


% remove duplicate

tmp_l = ifft(c_ls_FOV(:,:,1), NFFT, 1, 'symmetric');
tmp_l(NFFT_V+1:end,:) = 0;
tmp_l = fft(tmp_l,NFFT,1);
tmp_l = tmp_l(1:NFFT/2+1,:);
tmp_r = ifft(c_ls_FOV(:,:,2), NFFT, 1, 'symmetric');
tmp_r(NFFT_V+1:end,:) = 0;
tmp_r = fft(tmp_r,NFFT,1);
tmp_r = tmp_r(1:NFFT/2+1,:);
c_ls_FOV = cat(3,tmp_l,tmp_r);


tmp_l = ifft(c_mls_FOV(:,:,1), NFFT, 1, 'symmetric');
tmp_l(NFFT_V+1:end,:) = 0;
tmp_l = fft(tmp_l,NFFT,1);
tmp_l = tmp_l(1:NFFT/2+1,:);
tmp_r = ifft(c_mls_FOV(:,:,2), NFFT, 1, 'symmetric');
tmp_r(NFFT_V+1:end,:) = 0;
tmp_r = fft(tmp_r,NFFT,1);
tmp_r = tmp_r(1:NFFT/2+1,:);
c_mls_FOV = cat(3,tmp_l,tmp_r);

tmp_l = ifft(c_ls(:,:,1), NFFT, 1, 'symmetric');
tmp_l(NFFT_V+1:end,:) = 0;
tmp_l = fft(tmp_l,NFFT,1);
tmp_l = tmp_l(1:NFFT/2+1,:);
tmp_r = ifft(c_ls(:,:,2), NFFT, 1, 'symmetric');
tmp_r(NFFT_V+1:end,:) = 0;
tmp_r = fft(tmp_r,NFFT,1);
tmp_r = tmp_r(1:NFFT/2+1,:);
c_ls = cat(3,tmp_l,tmp_r);

tmp_l = ifft(c_mls(:,:,1), NFFT, 1, 'symmetric');
tmp_l(NFFT_V+1:end,:) = 0;
tmp_l = fft(tmp_l,NFFT,1);
tmp_l = tmp_l(1:NFFT/2+1,:);
tmp_r = ifft(c_mls(:,:,2), NFFT, 1, 'symmetric');
tmp_r(NFFT_V+1:end,:) = 0;
tmp_r = fft(tmp_r,NFFT,1);
tmp_r = tmp_r(1:NFFT/2+1,:);
c_mls = cat(3,tmp_l,tmp_r);





% Calculate for 'ls_fov' condition
[res_out.BRIR_ls_fov, res_out.BRTF_ls_fov, res_out.p_BSM_t_ls_fov] = ...
    calculateBRIR(c_ls_FOV, p_array_f, NFFT, s);

% Calculate for 'mls_fov' condition
[res_out.BRIR_mls_fov, res_out.BRTF_mls_fov, res_out.p_BSM_t_mls_fov] = ...
    calculateBRIR(c_mls_FOV, p_array_f, NFFT, s);

% Calculate for 'ls_normal' condition
[res_out.BRIR_ls_normal, res_out.BRTF_ls_normal, res_out.p_BSM_t_ls_normal] = ...
    calculateBRIR(c_ls, p_array_f, NFFT, s);

% Calculate for 'mls_normal' condition
[res_out.BRIR_mls_normal, res_out.BRTF_mls_normal, res_out.p_BSM_t_mls_normal] = ...
    calculateBRIR(c_mls, p_array_f, NFFT, s);

%fprintf('Finished all BRIR calculations.\n');

% Calculate reference HOA signal

pad = ceil(BSMobj.r_array(1)/BSMobj.c * fs * 100); % r/c*100 is just a bound on the length of the impulse response of the system from anm to p.
RIR_SH_t(end+1:end+1+pad, :) = 0; 

% FFT

%NFFT = 2^nextpow2(size(RIR_SH_t,1));
RIR_SH_f = fft(RIR_SH_t, NFFT, 1);

% remove negative frequencies
RIR_SH_f(NFFT/2+2:end,:)=[]; % perhaps a mistake
N_PW = sqrt(size(RIR_SH_f,2))-1;

tild_N_inf_mat = tildize(N_PW);
pad_center = (N_PW+1)^2;

anm_tilde_full = (RIR_SH_f*tild_N_inf_mat);

% load HRTF
hobj = sofaToEaro(convertStringsToChars(HRTF_path));
hobj.shutUp = true;
% hobj = hobj.toFreq(512);            % Transform from time to Freq domain
% hobj = hobj.toTime(512);            % Transform from time to Freq domain
% hobj.data = circshift(hobj.data,230,2);
hobj = hobj.resampleData(fs);
hobj = hobj.toSH(N_PW,'SRC');        % Transform from space to SH domain
hobj = hobj.toFreq(NFFT);            % Transform from time to Freq domain
hobj.data(:,NFFT/2+2:end,:) = [];    % get only the positive + DC freqncies
H_l_nm_full = hobj.data(:,:,1).';
H_r_nm_full = hobj.data(:,:,2).';


BRIR_REF_SH_L = sum(padarray(anm_tilde_full,[0 pad_center-size(anm_tilde_full,2)],'post').*padarray(H_l_nm_full,[0 pad_center-size(H_l_nm_full,2)],'post'),2);
BRIR_REF_SH_R = sum(padarray(anm_tilde_full,[0 pad_center-size(anm_tilde_full,2)],'post').*padarray(H_r_nm_full,[0 pad_center-size(H_r_nm_full,2)],'post'),2);

BRIR_REF_f = cat(2,BRIR_REF_SH_L,BRIR_REF_SH_R);



BRIR_REF_f(end+1:NFFT, :) = 0;

BRIR_REF_t = ifft(BRIR_REF_f, "symmetric");

BRIR_REF_t = BRIR_REF_t(1:size(res_out.BRIR_mls_normal,1),:);

%BRIR_REF_t = circshift(BRIR_REF_t,-1*(256 + 26),1);
BRIR_REF_t = circshift(BRIR_REF_t,NFFT_V/2 ,1);

%BRIR_REF_t = MatchEnergy(res_out.BRIR_mls_fov, BRIR_REF_t, fs, 50);

p_REF_t    = fftfilt(BRIR_REF_t, s);

% norm volume
%p_REF_t_norm = p_REF_t./max(max(abs(p_REF_t)));
p_REF_t_norm = p_REF_t;
%BRIR_REF_norm = BRIR_REF_t./max(max(abs(BRIR_REF_t)));
BRIR_REF_norm = BRIR_REF_t;


res_out.BRTF_ref    = BRIR_REF_f;
res_out.BRIR_ref    = BRIR_REF_norm;
res_out.p_BSM_t_ref = p_REF_t_norm;

% circshift to the front:
%res_out.BRIR_ref        = circshift(res_out.BRIR_ref,NFFT_V,1);
res_out.BRIR_mls_normal = circshift(res_out.BRIR_mls_normal,NFFT_V,1);
res_out.BRIR_mls_fov    = circshift(res_out.BRIR_mls_fov,NFFT_V,1);

%res_out.BRIR_mls_normal = AlignFilters_Time(res_out.BRIR_ref ,res_out.BRIR_mls_normal);
%res_out.BRIR_mls_fov    = AlignFilters_Time(res_out.BRIR_ref ,res_out.BRIR_mls_fov);




end

function sig_out = MatchEnergy(sig_ref, sig_in, fs, f_cutoff)
    % Matches energy based on LOW FREQUENCIES ONLY (< f_cutoff)
    % 1. Compute FFT
    nfft = 2^nextpow2(size(sig_ref,1));
    X_ref = fft(sig_ref, nfft);
    X_in  = fft(sig_in, nfft);
    
    % 2. Define Frequency Axis
    f_axis = linspace(0, fs, nfft)';
    
    % 3. Select Low Freq Bins
    idx_low = (f_axis <= f_cutoff) | (f_axis >= fs - f_cutoff);
    
    % 4. Compute Energy in Low Band
    E_ref_low = sum(abs(X_ref(idx_low, :)).^2, 1);
    E_in_low  = sum(abs(X_in(idx_low, :)).^2, 1);
    
    % 5. Calculate Gain Factor
    gains = mean(sqrt(E_ref_low ./ (E_in_low + 1e-12)));
    
    % 6. Apply Gain to Full Time Signal
    sig_out = sig_in .* gains;
end

function Perm=SHc2r(Nmax)

    % this code forms a permute matrix from the Normalized Complex Spherical Harmonics to
    % the Normalized Real Spherical Harmonics
    % Perm matrix hold the relation- Ynm_{Real} = Perm x Ynm_{Complex}

    Perm = zeros((Nmax+1)^2);
    sizeP = size(Perm,1);
    ind = 0;
    for n= 0:Nmax

        Perm((ind+1):(ind+1+(2*n+1)-1),(ind+1):(ind+1+(2*n+1)-1)) = miniSHc2r(n);
        ind = ind + (2*n +1);
    end

    Perm=conj(Perm);
   
end

function perm=miniSHc2r(n)

    % a help function for SHc2r, permuting for each given n.

    perm = zeros((2*n+1));
    sizeP = size(perm,1);
    perm((floor(sizeP/2)+1),(floor(sizeP/2)+1)) = 1;
    for ii= 1:(floor(sizeP/2))
        perm((floor(sizeP/2)+1+ii),(floor(sizeP/2)+1+ii)) = 1/sqrt(2)*(-1)^ii;%*(-1)^ii;
        perm((floor(sizeP/2)+1+ii),(floor(sizeP/2)+1-ii)) = 1/sqrt(2);
        perm((floor(sizeP/2)+1-ii),(floor(sizeP/2)+1-ii)) = -1/(sqrt(2)*1j);%*(-1)^ii;
        perm((floor(sizeP/2)+1-ii),(floor(sizeP/2)+1+ii)) = +1/(sqrt(2)*1j)*(-1)^ii;
    end
end


function [BRIR, BRTF, p_BSM_t_norm] = calculateBRIR(bsm_coeffs, p_array_f, NFFT, s)
% calculateBRIR computes the Binaural Room Impulse Response (BRIR) and 
% Transfer Function (BRTF) from Binaural Signal Matching (BSM) coefficients.
%
% This function correctly constructs a conjugate-symmetric spectrum before
% the IFFT to produce a real-valued, causal impulse response, avoiding the
% "mirror image" artifact.
%
% Inputs:
%   bsm_coeffs  - 3D array of BSM coefficients [nFreq, nMicrophones, 2 ears].
%   p_array_f   - 2D array of the frequency-domain microphone array signals 
%                 [nFreq, nMicrophones].
%   NFFT        - The desired FFT size for the impulse response.
%   s           - The source signal to be convolved with the BRIR using fftfilt.
%
% Outputs:
%   BRIR        - The calculated Binaural Room Impulse Response [NFFT, 2].
%   BRTF        - The calculated Binaural Room Transfer Function (positive 
%                 frequencies only) [nFreq, 2].
%   p_BSM_t_norm- The normalized time-domain signal after convolving the
%                 source 's' with the calculated BRIR.

    % Initialize arrays to hold the left and right channel results.
    BRIR_channels = zeros(NFFT, 2);
    BRTF_channels = zeros(size(p_array_f, 1), 2);

    % Process both ears (channels). The BSM coefficients should have the
    % 3rd dimension of size 2.
    for q = 1:2
        
        % --- Core Calculation with Hermitian Symmetry Fix ---
        
        % 1. Calculate the BRTF for the positive frequencies (from DC to Nyquist).
        % This is the transfer function for one ear.
        BRTF_pos = sum(conj(bsm_coeffs(:,:,q)) .* p_array_f, 2);
        %BRTF_pos = sum(bsm_coeffs(:,:,q) .* p_array_f, 2);
        
        % Store this positive-frequency transfer function for the final output.
        BRTF_channels(:, q) = BRTF_pos;
        
        % 2. Get the number of positive frequency points.
        nFreq = size(BRTF_pos, 1);
        
        % 3. Construct the negative frequency components to ensure the full 
        % spectrum has conjugate (Hermitian) symmetry. This is CRITICAL to 
        % get a real-valued signal after the IFFT.
        if mod(NFFT, 2) == 0 
            % Case for an EVEN NFFT length.
            % The unique points are [DC, f1, ..., f_nyquist]. Size is NFFT/2 + 1.
            if nFreq ~= (NFFT/2 + 1)
                warning('calculateBRIR:NFFT_mismatch', ...
                        'For even NFFT, expected %d freq points, but got %d.', (NFFT/2 + 1), nFreq);
            end
            % The negative frequencies are the complex conjugate of the positive ones,
            % in reverse order. We exclude the DC (index 1) and Nyquist (index end)
            % components as they do not have a conjugate pair.
            BRTF_neg = conj(BRTF_pos(end-1:-1:2, :));
        else 
            % Case for an ODD NFFT length.
            % The unique points are [DC, f1, ..., f_max]. Size is (NFFT+1)/2.
             if nFreq ~= ((NFFT+1)/2)
                warning('calculateBRIR:NFFT_mismatch', ...
                        'For odd NFFT, expected %d freq points, but got %d.', ((NFFT+1)/2), nFreq);
            end
            % For odd NFFT, there is no Nyquist frequency. We only exclude DC.
            BRTF_neg = conj(BRTF_pos(end:-1:2, :));
        end
        
        % 4. Assemble the full, conjugate-symmetric spectrum.
        full_BRTF = [BRTF_pos; BRTF_neg];
        
        % 5. Perform the IFFT. The 'symmetric' flag tells ifft that the
        % input is conjugate-symmetric, ensuring the output is purely real.
        BRIR_channels(:, q) = ifft(full_BRTF, NFFT, 1, 'symmetric');
    end

    % Assign the final outputs from the processed left and right channels.
    BRIR = BRIR_channels;
    
    % %BRIR(1:nFreq,:) = 0;
    % tmp = BRIR(:,1);
    % BRIR(:,1) = BRIR(:,2);
    % BRIR(:,2) = tmp;

    BRTF = BRTF_channels;

    % Convolve the resulting BRIR with the source signal.
    p_BSM_t = fftfilt(BRIR, s);
    
    % Normalize the convolved signal to a peak amplitude of 1.
    % This avoids division by zero if the signal is silent.
    max_abs_val = max(abs(p_BSM_t), [], 'all');
    if max_abs_val > 0
        p_BSM_t_norm = p_BSM_t ./ max_abs_val;
    else
        p_BSM_t_norm = p_BSM_t; % Signal is zero, return as is.
    end
end


function c_interp = interpolateBSM(c_orig, f_orig, f_interp)
% interpolateBSM robustly interpolates BSM coefficients by operating on
% magnitude and unwrapped phase separately to avoid time-domain aliasing.
%
% Inputs:
%   c_orig      - 3D array of the original, uninterpolated BSM coefficients
%                 [nFreq_orig, nMicrophones, 2 ears].
%   f_orig      - The original frequency vector corresponding to c_orig.
%   f_interp    - The target frequency vector for interpolation.
%
% Output:
%   c_interp    - The interpolated BSM coefficients [nFreq_interp, nMicrophones, 2].

    [~, I, nEars] = size(c_orig);
    nFreq_interp = length(f_interp);
    c_interp = zeros(nFreq_interp, I, nEars, 'like', 1j); % Ensure complex output

    for i = 1:I
        for q = 1:nEars
            % Extract one vector of complex coefficients
            coeffs_vec = squeeze(c_orig(:, i, q));
            
            % --- Convert to Magnitude and Unwrapped Phase ---
            % This is the key step to allow for smooth interpolation.
            mag = abs(coeffs_vec);
            % unwrap() removes phase jumps (e.g., from +pi to -pi)
            phase = unwrap(angle(coeffs_vec));
            
            % --- Interpolate Magnitude and Phase Separately ---
            interp_mag = interp1(f_orig, mag, f_interp, 'pchip', 'extrap');
            interp_phase = interp1(f_orig, phase, f_interp, 'pchip', 'extrap');
            
            % --- Reconstruct the Complex Coefficients ---
            c_interp(:, i, q) = interp_mag .* exp(1j * interp_phase);
        end
    end
end


function c_new = convertBSM(c_orig, NFFT)
% convertBSM converts BSM coefficients from a circular convolution
% representation to a linear one to prevent time-domain aliasing.
%
% This is achieved by transforming the coefficients to the time-domain,
% zero-padding them to the new desired length, and transforming back to the
% frequency domain.
%
% Inputs:
%   c_orig - Original BSM coefficients [L_orig, nMics, 2 ears].
%   NFFT   - The target FFT length for the linear convolution.
%
% Output:
%   c_new  - New BSM coefficients correctly representing linear convolution
%            for the given NFFT size [NFFT/2 + 1, nMics, 2 ears].

    % Get original length from the input coefficients.
    L_orig = size(c_orig, 1);
    
    % Determine the number of positive frequency points for the new FFT size.
    nFreq_new = NFFT/2 + 1;
    
    % 1. IFFT the original coefficients to get the aliased time-domain IR.
    ir_aliased = ifft(c_orig, L_orig, 1, 'symmetric');
    
    % 2. Create a new, larger array for the zero-padded time-domain IR.
    %    Initialize with zeros to ensure proper padding.
    ir_padded = zeros([NFFT, size(ir_aliased, 2), size(ir_aliased, 3)]);
    
    %    Copy the original IR into the start of the new padded array.
    ir_padded(1:L_orig, :, :) = ir_aliased;
    
    % 3. FFT the padded IR to get the new frequency-domain coefficients.
    c_fft = fft(ir_padded, NFFT, 1);
    
    % 4. Return only the positive frequency components, which is what the
    %    'symmetric' flag requires for ifft.
    c_new = c_fft(1:nFreq_new, :, :);
end