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
% B/S/A BEMA Spatial Anti Aliasing R13-0603
%
%       BEMA Bandwidth Extension for Microphone Arrays
%       -> AES Convention 2012, San Francisco USA
%          Convention Paper 8751
%
% Copyright (C)2012-2013 Benjamin Bernschütz
%
%
% This file is part of the SOFiA toolbox under GNU General Public License
%
%
% Pnm = sofia_bsa(Pnm, ctSig, dn, transition, avgBandwidth, fade)
% ------------------------------------------------------------------------
% Pnm               Alias-Free Spatial Fourier Coefficients
%
% ------------------------------------------------------------------------
% Pnm               Spatial Fourier Coefficients
% ctSig             Signal of the center microphone
% dn                Radial filters for the current array configuration
% transition        Highest stable bin
%                   Approx: transition = (NFFT/FS+1) * (N*c)/(2*pi*r)
% avgBandwidth      Averaging Bandwidth in oct
% fade              Fade over {true} or hard cut {false}


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


function Pnm = sofia_bsa(Pnm, ctSig, dn, transition, avgBandwidth, fade)

disp('SOFiA B/S/A - BEMA Spatial Anti Aliasing R13-0306');

startaverage = floor(transition/(2^(avgBandwidth))); %Averaging BW
dn_inverted  = 1./dn;                                %Inverted Radial Filters

SpatialImage = zeros(size(Pnm,1),1);
AveragePower          = 0;
avgBins               = length(startaverage:transition);
N                     = sqrt(size(Pnm,1))-1;

cnt=0;

%Extraction of the Spatial Image
for n=0:N
    for m=-n:n
        cnt=cnt+1;
        for bincounter = startaverage:transition
            compensatedPnmBin  = Pnm(cnt,bincounter) / dn_inverted(n+1,bincounter);
            SpatialImage(cnt)  = SpatialImage(cnt) + compensatedPnmBin;
            AveragePower       = AveragePower      + abs(compensatedPnmBin).^2;
        end
    end
end

%Normalization of the Image
SumEnergyOfAllModesInFP  = sum(abs(SpatialImage).^2);
powerNormalizationFactor = AveragePower/avgBins;
SpatialImage    = SpatialImage*sqrt(powerNormalizationFactor/SumEnergyOfAllModesInFP);

%Normalization of the Center Signal
ctSigAverage = mean(abs(ctSig(startaverage:transition)).^2);
ctSig        = ctSig/sqrt(ctSigAverage);

%Add spectral information from Center Signal to Spatial Image
Pnm_bema = zeros(size(Pnm));
cnt = 0;

for n=0:N
    for m=-n:n
        cnt=cnt+1;
        for f=startaverage:size(Pnm,2)
            Pnm_bema(cnt, f) = SpatialImage(cnt) * dn_inverted(n+1,f)* ctSig(f);
        end
    end
end

%Phase correction
phaseoffset = angle(Pnm(1, transition)/Pnm_bema(1,transition));
Pnm_bema = Pnm_bema.*exp(1i*phaseoffset);

%Replace high bins with BEMA-coefficients
Pnm = [Pnm(:, 1:transition-1) Pnm_bema(:, transition:end)];

%Fade over original coefficients and BEMA-coefficients
if fade
    PnmFade_O = Pnm(:, startaverage:transition-1);
    PnmFade_S = Pnm_bema(:, startaverage:transition-1);
    
    fader    = linspace(0,1,size(PnmFade_O,2));
    fadeup   = repmat(fader,size(PnmFade_O,1),1); %Make Matrix UP
    fadedown = fliplr(fadeup);                    %Make Matrix DOWN
    
    PnmFade_O = PnmFade_O.*fadedown;
    PnmFade_S = PnmFade_S.*fadeup;
    
    PnmFade = PnmFade_O + PnmFade_S;
    Pnm(:, startaverage:transition-1) = PnmFade;
end

