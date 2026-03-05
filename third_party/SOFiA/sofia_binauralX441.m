% /// ASAR-MARA Research Group
%
% Cologne University of Applied Sciences
% Berlin University of Technology
% University of Rostock
% Deutsche Telekom Laboratories
% WDR Westdeutscher Rundfunk
% IOSONO GmbH
% 
% SOFiA sound field analysis
% 
% BinauralX441 - Binaural Synthesis R13-0306
% 
% Copyright (C)2011-2013 Benjamin Bernschütz   
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
%
% [BinauralL, BinauralR] = sofia_binauralX441(Pnm, dn, headRot, [cGridType], ...
%                          [hpfOn], [hpTaps], [ncGap])
% ------------------------------------------------------------------------     
% BinauralL/R        Binaural Frequency domain FFT data for multiple ch.
%                    Columns : Index / Channel
%                    Rows    : FFT data (frequency domain)
% ------------------------------------------------------------------------              
% Pnm                Spatial Fourier Coefficients
% 
% dn                 Modal Array Filters from SOFiA M/F
%                   
% headRot            Head rotation(s) in RAD
%                  
% [cGridType]        Composite grid type 
%                    0 Gauss-Legendre quadrature (#default)
%                    1 Lebedev quadrature 
%                    
% [hpfOn]            30Hz highpass filter 
%                    true  : HPF on (#default)
%                    false : HPF off
%
% [hpTaps]           Taps for the highpass filter, 2048 (#default)
%
% [ncGap]            Non-causality FIR gap (radial filters), 16 (#default)
%                     -> Entails latency of ncGap samples
%
%
% The function synhtesizes binaural signals from spatial fourier coeff- 
% icients Pnm. The binaural cues are generated using the HRTFs of a Neumann 
% KU100 artificial head. The function combines plane wave decomposition on 
% a matched order composite grid and spatial downsampling of high order  
% (N=35) interpolated HRTFs. The binaural signals can be directly calculated 
% for different azimutal head rotations in a single run (e.g. to generate 
% full dynamic binaural synthesis datasets).
%
% ! IMPORTANT ADVICE:
%
%  I.  Only signals with a sampling rate of 44100Hz supported (HRTFs)
%  II. Initial NFFT must be at least twice the length of the IR
%    
%
% Dependencies: - sofia_pdc               
%               - sofia_itc
%               - sofia_lebedev or sofia_gauss
%               - SOFiA KU100 HRIR Set stored in: "sofia_HRIR.mat"
%               - MATLAB Signal Processing Toolbox for the HPF   
%

% CONTACT AND LICENSE INFORMATION:
%
% /// ASAR-MARA Research Group
%
%     [1] Cologne University of Applied Sciences
%     [2] Berlin University of Technology
%     [3] Deutsche Telekom Laboratories
%     [4] WDR Westdeutscher Rundfunk
%     [5] University of Rostock
%     [6] IOSONO GmbH
%
% SOFiA sound field analysis
%
% Copyright (C)2011-2013 Benjamin Bernschütz [1,2] et al.(§)
%
% Contact -------------------------------------
% Cologne University of Applied Sciences
% Institute of Communication Systems
% Betzdorfer Street 2
% D-50679 Germany (Europe)
%
% phone +49 221 8275 -2496
% mail  benjamin.bernschuetz@fh-koeln.de
% ---------------------------------------------
%
% This file is part of the SOFiA sound field analysis toolbox
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
% (§) Christoph Pörschmann [1]   christoph.poerschmann 'at' fh-koeln.de
%     Stefan Weinzierl     [2]   stefan.weinzierl 'at' tu-berlin.de
%     Sascha Spors         [5]   sascha.spors 'at' uni-rostock.de


function [BinauralL, BinauralR] = sofia_binauralX441(Pnm, dn, headRot, ...
                                  cGridType, hpfOn, hpTaps, ncGap)

if nargin < 7
    ncGap = 16; 
end

if nargin < 6
    hpTaps = 2048; 
end

if nargin < 5
    hpfOn = true; 
end

if nargin < 4
    cGridType = 0;
end

clc
fprintf('SOFiA BinauralX441 - Binaural Synthesis R13-0306\n\n');
fprintf(' -> Only signals with a sampling rate of 44100Hz supported!\n')
fprintf(' -> Initial NFFT must be at least twice the length of the IR!\n\n')
fprintf('Be patient - this may take a while...\n\n')
pause(1.5)

if ~exist('sofia_HRIR441.mat','file')
    error('Binaural HRIR dataset not accessible! File: sofia_HRIR441.mat');
else
    load sofia_HRIR441
end

[compositeGrid, Npwd] = getCompositeGrid(Pnm, cGridType);

Y       = sofia_pdc(Npwd, compositeGrid(:,1:2), Pnm, dn);
weights = repmat(compositeGrid(:,3),1,size(Y,2));
Y       =  Y .* weights * (Npwd+1)^2;

[BinauralL, BinauralR] = binauralize(Hl_nm, Hr_nm, Y, compositeGrid, headRot, hpfOn, hpTaps, ncGap);



function [compositeGrid, Npwd] = getCompositeGrid(Pnm, cGridType)

Npwd = sqrt(size(Pnm,1))-1;

if cGridType
    
    lebNodes = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, ...
        350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, ...
        3074, 3470, 3890, 4334, 4802, 5294, 5810];
    
    lebOrders  = floor(sqrt(lebNodes/1.3)-1);
    matchOrder = abs(lebOrders-Npwd);
    orderIndex = find(matchOrder==min(matchOrder));
    
    if ~matchOrder(orderIndex)
        fprintf(['>>> composite grid: Lebedev, ', num2str(lebNodes(orderIndex)),' Nodes\n\n']);
        compositeGrid = sofia_lebedev(lebNodes(orderIndex),0);
        return
    else
        warning('No matching Lebedev composit grid available: Switched to Gauss grid.\n');
    end
end

fprintf(['>>> composite grid: Gauss, ', num2str(2*(Npwd+1)^2),' Nodes\n\n']);
compositeGrid = sofia_gauss(2*(Npwd+1), (Npwd+1),0);



function [BinauralL, BinauralR] = binauralize(Hl_nm, Hr_nm, Y, compositeGrid, headRot, hpfOn, hpTaps, ncGap)

NFFT = (2*(size(Y,2)-1));

Y = conj([Y(:,:), conj(fliplr(Y(:,2:end-1)))])';
y = real(ifft(Y));
y = circshift(y,ncGap); % cyclic move for ncGap

u      = 0:2*ncGap-1;   % window full block (radial filters)
winFkt = 0.5+0.5*cos(2*pi*(u-((2*ncGap-1)/2))/(2*ncGap-1));
winFkt = [winFkt(1:end/2), ones(1,NFFT/2), winFkt(end/2+1:end), zeros(1,NFFT/2-2*ncGap)]';
winFkt = repmat(winFkt,1,size(y,2));
y      = y.*winFkt;

if hpfOn %Recalc required NFFT
    NFFT = NFFT/2 + 2*ncGap + (size(Hl_nm,2)-1) + hpTaps; 
else
    NFFT = NFFT/2 + 2*ncGap + (size(Hl_nm,2)-1); %reserve 128taps only for HRTF 
end

if size(y,1) > NFFT
    y = y(1:NFFT,:);
end

Y = fft(y,NFFT);

BinauralL = zeros(length(headRot),NFFT);
BinauralR = zeros(length(headRot),NFFT);

fprintf('\n');

for hrPointer = 1:length(headRot)
    
    fprintf(['Head orientation ', num2str(mod(headRot(hrPointer)*180/pi,360)),'°\n']);
    hrtfRotAZ           = mod(compositeGrid(:,1)-headRot(hrPointer),2*pi);
    Hl                  = sofia_itc(Hl_nm, [hrtfRotAZ compositeGrid(:,2)]);
    Hr                  = sofia_itc(Hr_nm, [hrtfRotAZ compositeGrid(:,2)]);
    
    Hl   = conj([Hl(:,:), conj(fliplr(Hl(:,2:end-1)))])';
    Hr   = conj([Hr(:,:), conj(fliplr(Hr(:,2:end-1)))])';
    hl   = real(ifft(Hl));
    hr   = real(ifft(Hr));
    Hl   = fft(hl,NFFT);
    Hr   = fft(hr,NFFT);

    BinL = Y.*Hl;
    BinL = sum(BinL,2);
    BinauralL(hrPointer, :) = BinL;
    BinR = Y.*Hr;
    BinR = sum(BinR,2);
    BinauralR(hrPointer, :) = BinR;    
    
end

if hpfOn
    
    hpf       = 30;    
    fs        = 44100;
    hpf       = fir1(hpTaps,hpf/(fs/2),'high');
    HPF       = fft(hpf, 2*hpTaps);
    img       = imag(hilbert(log(abs(HPF))));
    hpf       = real(ifft(abs(HPF) .* exp(-1i*img)));
    hpf       = hpf(1:end/2);
    winLength = round(hpTaps/2);
    c         = 0:(round(winLength))-1;
    winFkt    = 0.5+0.5*cos(2*pi*(c-((winLength-1)/2))/(winLength-1));
    winFkt    = winFkt(end/2+1:end);
    winFkt    = [ones(1,hpTaps-size(winFkt,2)), winFkt];
    hpf       = hpf.*winFkt;
    HPF       = fft(hpf, NFFT);
    HPF       = repmat(HPF,size(BinauralL,1),1);
    BinauralL = BinauralL.*HPF;
    BinauralR = BinauralR.*HPF;
    
end

BinauralL = BinauralL(:,1:end/2+1);
BinauralR = BinauralR(:,1:end/2+1);


