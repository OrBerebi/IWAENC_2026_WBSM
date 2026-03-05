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
% CAD Reflection Analysis: Equal Distance Matrix 
%                          Generator Wrapper R13-01-24
%
% Copyright (C)2013 by Maximilian Rühl and Benjamin Bernschütz
%
% This file is part of the SOFiA toolbox under GNU General Public License
%
% eqDistMTX = makeEqDistMTX(timeData, ac, frequency, timeLimit, lebOrder,
%                           blockSize, foreShift, N, maxAmp, fftOversize)
% ------------------------------------------------------------------------
% eqDistMTX         Struct with following fields:
%                   * .mtxData
%                   * .ac
%                   * .frequency
%                   * .timeLimit
%                   * .lebOrder
%                   * .blockSize
%                   * .foreshift
%                   * .N
%                   * .maxAmp
%                   * .fftOversize
% ------------------------------------------------------------------------
% timeData          Struct with minimum fields:
%                   * .impulseResponses     [Channels X Samples]
%                   * .FS
%                   * .radius               Array radius
%                   * .averageAirTemp       Temperature in [°C]
%                   (* .centerIR            [1 x Samples] )
%
% ac                Array Configuration:
%                   0  Open sphere with pressure transducers (NO plc!)
%                   1  Open sphere with cardioid transducers
%                   2  Rigid sphere with pressure transducers
%                   3  Rigid sphere with cardioid transducers
%                   4  Dual open sphere with pressure transducers
%
% frequency         Frequency in Hz, [Default 2000Hz]
%
% timelimit         Time Limit in Seconds, [Default: 100ms]
%
% lebOrder          Lebedev Grid Order (Number of Angles) [Default 3074]
%
% blockSize         FFT Blocksize [Default 64]
%
% foreShift         Window ForeShift [Default blockSize/2]
%
% N                 Spherical Decomposition Order [Default 5]
%
% maxAmp            Maximum modal amplification limit in dB
%
% FFToversize       FFT Oversize [Default 1]
%
% This Function generates a sound field matrix based on a Lebedev
% decomposition grid.


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
%     Maximilian Rühl      [1]   maxruehl 'at msn.com


function eqDistMTX = sofia_makeEqDistMTX(timeData, ac, frequency, timeLimit, lebOrder, blockSize, foreShift, N, maxAmp, fftOversize)

%CONFIGURATION
if nargin < 2,  ac           = 2;           end % Array configuration (compare SOFiA M/F)
if nargin < 3,  frequency    = 2000;        end % Frequency in Hz, [Default 2000Hz]
if nargin < 4,  timeLimit    = 0.100;       end % Time Limit in Seconds, [Default: 100ms]
if nargin < 5,  lebOrder     = 3074;        end % Lebedev Grid Order (Number of Angles) [Default 3074]
if nargin < 6,  blockSize    = 64;          end % FFT Blocksize [Default 64]
if nargin < 7,  foreShift    = blockSize/2; end % Window ForeShift [Default blockSize/2]
if nargin < 8,  N            = 5;           end % Spherical Decomposition Order [Default 5]
if nargin < 9,  maxAmp       = 10;          end % Maximum modal amplification in dB
if nargin < 10, fftOversize  = 1;           end % FFT Oversize [Default 1]

disp('SOFiA Equal Distance Matrix Generator Wrapper R13-01-24');

lengthOfIRs       = size(timeData.impulseResponses,2);
lastPossibleBlock = lengthOfIRs - blockSize;
timeLimitSamples  = timeLimit*timeData.FS;

if timeLimitSamples<lastPossibleBlock;
    lastPossibleBlock = timeLimitSamples;
end

frames = floor(lastPossibleBlock / foreShift);

frameResolution = foreShift / timeData.FS;

frameVector = zeros(frames,1);
timeVector  = zeros(frames,1);

for i=1:frames
    frameVector(i) = i*foreShift;
    timeVector(i)  = i*frameResolution;
end

% Generate first Vector and Radial Filters
% Calculation of kr only (fftData is not used):

[notUsed, kr, f] = sofia_fdt(timeData, fftOversize, frameVector(1), frameVector(1)+blockSize);
dn               = sofia_mf(N, kr, ac, maxAmp);

krIndex = round(frequency/((timeData.FS/2)/(blockSize*fftOversize)))+1;

disp(['Frequency       : ',num2str(round(10*f(krIndex))/10),'Hz']);
disp(['Time Resolution : ',num2str(frameResolution*1000),'ms']);

quad = sofia_lebedev(lebOrder,0);
quadPDC = quad(:,1:2);

for i = 1:size(frameVector,1)
    
    pause(0.1)
    fftData  = sofia_fdt(timeData, fftOversize, frameVector(i), frameVector(i)+blockSize);
    Pnm      = sofia_stc(N, fftData, timeData.quadratureGrid);
    
    Y = sofia_pdc(N, quadPDC, Pnm(:,krIndex), dn(:,krIndex));
    Y = abs(Y);
    
    mtxData(:,:,i) = [Y quad];   
    
end

eqDistMTX.mtxData      = mtxData;
eqDistMTX.ac           = ac;
eqDistMTX.frequency    = frequency;
eqDistMTX.timeLimit    = timeLimit;
eqDistMTX.lebOrder     = lebOrder;
eqDistMTX.blockSize    = blockSize;
eqDistMTX.foreShift    = foreShift;
eqDistMTX.N            = N;
eqDistMTX.maxAmp       = maxAmp;
eqDistMTX.fftOversize  = fftOversize;








