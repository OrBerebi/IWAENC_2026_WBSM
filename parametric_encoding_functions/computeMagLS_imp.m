function [H_l_nm_MagLS, H_r_nm_MagLS] = computeMagLS_imp(H, f_vec, N, Y_N, cutOffFreq, Y_high,IACC_correction)
% function computeMagLS computes the MagLS representation of the HRTFs
%
% Implementation based on Ambisonics book (https://link.springer.com/book/10.1007/978-3-030-17207-7)
%   Eqs.(4.57)-(4.59), for N<3 until (4.63)
%
% Inputs:
%       H          : HRTF in space/frequncy domain [# of directions X # of frequency bins X 2 (left/right)]
%       f_vec      : Frequencies vector [Hz]
%       N          : SH order
%       Y_N        : Complex normalized SH functions matrix of order N [(N+1)^2 X # of directions ]
%       cutOffFreq : Frequency cut-off for MagLS [Hz] (set at 2KHz by deafult)
%       Y_high     : Complex normalized SH functions matrix of high order M>>N (used for Diffuse-Field Covariance Constraint, when N<3)
% Output:
%       H_l_nm_MagLS : SH coefficient of left HRTF [# of frequency bins X (N+1)^2]
%       H_r_nm_MagLS : SH coefficient of right HRTF [# of frequency bins X (N+1)^2]
%
% Written by Zamir Ben-Hur
% Nov 2020
%
% Modified by Or B September 2021

% initials
if nargin < 6 % don't use Diffuse-Field Covariance Constraint
    Y_high = [];
    IACC_correction = true;
end

% calc start and end transition freq band
k = 0.5; % the length [in octave] of the transition (half)
%fc_low = cutOffFreq*2^(-k/2);
fc_low = cutOffFreq;
fc_high = cutOffFreq*2^(k);

% find cut-off freq index
srcharray = abs(f_vec-fc_low);
cutOffFreqInd_low = find(srcharray==min(srcharray));
srcharray = abs(f_vec-fc_high);
cutOffFreqInd_high = find(srcharray==min(srcharray));

% pseudo-inverse of matrix Y
%degrees = size(Y_N,2);
%[gridData_Inf, ~, ~]    = sofia_lebedev(degrees,0);
%ph = gridData_Inf(:,1); 
%th = gridData_Inf(:,2);
%a  = gridData_Inf(:,3)*4*pi;
%pY = diag(a)*Y_N';

pY = pinv(Y_N); 

H_l = squeeze(double(H(:,:,1))).'; % left HRTF
H_r = squeeze(double(H(:,:,2))).'; % right HRTF

H_l_nm_MagLS = zeros(size(H,2),(N+1)^2);
H_r_nm_MagLS = zeros(size(H,2),(N+1)^2);

% Compute MagLS

% for freq < cutOffFreq, standard SH transform (Eq. 4.57)
H_l_nm_MagLS(1:cutOffFreqInd_low-1,:) = H_l(1:cutOffFreqInd_low-1,:)*pY; % SH transform - left ear
H_r_nm_MagLS(1:cutOffFreqInd_low-1,:) = H_r(1:cutOffFreqInd_low-1,:)*pY; % SH transform - right ear



% for fc_low < freq < fc_high, MagLS (Eqs. 4.58, 4.59)
for freqInd = cutOffFreqInd_low : cutOffFreqInd_high
    %alpha = 0.5 + k/log(2) *  log(f_vec(freqInd)/cutOffFreq); % smoothing weight
    alpha = (log2(f_vec(freqInd)) - log2(fc_low)) / (log2(fc_high) - log2(fc_low));
    % left ear
    phi_est_l = angle(Y_N.'*H_l_nm_MagLS(freqInd-1,:).'); % Eq. (4.58)
    H_l_nm_MagLS(freqInd,:) = alpha*(abs(H_l(freqInd,:)) .* exp(1i*phi_est_l.')) * pY + (1- alpha)*(H_l(freqInd,:)*pY); %  Eq. (4.59)
    
    % right ear
    phi_est_r = angle(Y_N.'*H_r_nm_MagLS(freqInd-1,:).'); % Eq. (4.58)
    H_r_nm_MagLS(freqInd,:) = alpha*(abs(H_r(freqInd,:)) .* exp(1i*phi_est_r.')) * pY + (1- alpha)*(H_r(freqInd,:)*pY); %  Eq. (4.59)
end

%figure
% for freq >= cutOffFreq, MagLS (Eqs. 4.58, 4.59)
%tmp = abs(H_l(cutOffFreqInd_high,:));
for freqInd = cutOffFreqInd_high+1 : size(H,2)
    % left ear
    phi_est_l = angle(Y_N.'*H_l_nm_MagLS(freqInd-1,:).'); % Eq. (4.58)
    %plot(phi_est_l)
    %drawnow
    %tmp = ones(size(H_l(freqInd,:)));
    H_l_nm_MagLS(freqInd,:) = (abs(H_l(freqInd,:)) .* exp(1i*phi_est_l.')) * pY; %  Eq. (4.59)
    %H_l_nm_MagLS(freqInd,:) = (tmp .* exp(1i*phi_est_l.')) * pY; %  Eq. (4.59)
    % right ear
    phi_est_r = angle(Y_N.'*H_r_nm_MagLS(freqInd-1,:).'); % Eq. (4.58)
    H_r_nm_MagLS(freqInd,:) = (abs(H_r(freqInd,:)) .* exp(1i*phi_est_r.')) * pY; %  Eq. (4.59)
end

%for N<3 use Diffuse-Field Covariance Constraint
if IACC_correction
if N <= 3 && ~isempty(Y_high)
    [H_l_nm_MagLS,H_r_nm_MagLS] = diff_field_eq(Y_high,H,H_l_nm_MagLS,H_r_nm_MagLS);

end
end
end