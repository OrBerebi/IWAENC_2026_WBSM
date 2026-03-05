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
% A/W/G/N Additive White Gaussian Noise Generator R13-0306
% 
% Adds White Gaussian Noise of approx. 16dB crest to an fftData block in
% order to simulate equipment noise (from transducers, amplifiers,...).  
%
% Copyright (C)2013 Benjamin Bernschütz
%                   rockzentrale 'AT' me.com
%                   +49 171 4176069 Germany
%
% This file is part of the SOFiA toolbox under GNU General Public License
%
%
% noisyFftData = sofia_awgn(fftData, noiseLevel)
% -----------------------------------------------------------------------
%
% noisyFftData  Output fftData block including white gaussian noise 
%
% fftData       Input fftData block (e.g. from F/D/T or S/W/G)
% noiseLevel    Average noise Level in dB (#Default: -80dB)


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


function noisyFftData = sofia_awgn(fftData, noiseLevel)

if nargin < 2
   noiseLevel = -80;
end

disp('SOFiA A/W/G/N - Additive White Gaussian Noise Generator R13-0306');

dimFactor       = 10^(noiseLevel/20);
NFFT            = size(fftData,2)*2-2;
nNoise          = randn(size(fftData,1),NFFT);
nNoise          = nNoise/mean(mean(abs(nNoise)));
nNoise          = dimFactor*nNoise;
nNoiseSpectrum  = fft(nNoise,[],2);
nNoiseSpectrum  = nNoiseSpectrum(:,1:end/2+1);
noisyFftData    = fftData + nNoiseSpectrum;