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
% CAD Reflection Analysis: PK/D Peak Detection R13-01-24
%
% Copyright (C)2013 by Maximilian Rühl
%
% This file is part of the SOFiA toolbox under GNU General Public License
%
% detPeaks = sofia_peakDetection(eqDistMTX, dBsensitivity, surfaceRange,
%                                                   timeRange, dBthreshold)
% ------------------------------------------------------------------------
% detPeaks          All detected maxima within the eQDistMTX.mtxData
%                   Listing of local maxima where each row represents a
%                   local maximum containing the following properties:
%                   1. Respective time slice number
%                   2. Time in ms (Refers to the center of the FFT block)
%                   3. Azimut in DEG    (DEG: Orientation Purpose)
%                   4. Elevation in DEG (DEG: Orientation Purpose)
%                   5. Intensity        (Ref. to the absolute maximum)
%                   6. Intensity in dB  (Ref. to the absolute maximum)
%                   7. Azimut in RAD
%                   8. Elevation in RAD
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
%
% dBsensitivity     Threshold value referring to the mean intensity level
%                   within a time slice. A maximum candidate has to exceed
%                   this threshold to be valid.
%
% surfaceRange      Number of nodes of the Lebedev grid that are compared
%                   to a maximum candidate.
%
% timeRange         Number of time slices (+/-) that are compared to the
%                   maximum candidate.
%
% dBthreshold       Absolute intensity threshold. A maximum candidate has
%                   to exceed this threshold to be valid.
%
% This Function detects local maxima within a sound field matrix generated
% with makeEqDistMTX.



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
%     Maximilian Rühl      [1]   maxruehl 'at' msn.com


function detPeaks = sofia_pkd(eqDistMTX, dBsensitivity, surfaceRange, timeRange, dBthreshold)

if nargin < 5, dBthreshold   = 60; end
if nargin < 4, timeRange     = 1;  end
if nargin < 3, surfaceRange  = 20; end
if nargin < 2, dBsensitivity = 12; end

if nargin < 1
    error('eqDistMTX needed');
end

disp('SOFiA PK/D Peak Detection R13-03-06');

[Xg,Yg,Zg] = sph2cart(eqDistMTX.mtxData(:,2,1),eqDistMTX.mtxData(:,3,1)-pi/2,1);
lebSize    = size(eqDistMTX.mtxData,1);
timeSize   = size(eqDistMTX.mtxData,3);

distances  = zeros(lebSize,2,lebSize);
detPeaks  = [0,0,0,0,0,0,0,0]; % Init row

absoluteMaximum = max(max(eqDistMTX.mtxData(:,1,:)));

threshold   = absoluteMaximum*(10^(-dBthreshold/20));
sensitivity = 10^(dBsensitivity/20);

timeEnergy = mean(eqDistMTX.mtxData(:,1,:));

%Distances between the nodes of the Lebedev-Grid
for i=1:lebSize
    distances(:,1,i) = sqrt((Xg-Xg(i)).^2+(Yg-Yg(i)).^2+(Zg-Zg(i)).^2);
    distances(:,2,i) = 1:lebSize;
end

%Sorting the Distances
for i=1:lebSize
    distances(:,:,i) = sortrows(distances(:,:,i),1);
end

%progresss bar
progress = 1;
step     = timeSize/40;

for time = 1+timeRange:timeSize-timeRange
    
    % progress bar
    if time >= step*progress
        fprintf('|');
        progress=progress+1;
    end
    
    for i = 1:lebSize
        maximumFlag = 1;
        for compareTime = time-timeRange:time+timeRange
            if maximumFlag == 0, break, end
            for comparei = 1:surfaceRange
                if eqDistMTX.mtxData(i,1,time) < eqDistMTX.mtxData(distances(comparei,2,i),1,compareTime) %Comparison between 2 intensities
                    maximumFlag = 0;
                    break
                end
            end
        end
        if maximumFlag == 1
            if eqDistMTX.mtxData(i,1,time) >= sensitivity*timeEnergy(time) && eqDistMTX.mtxData(i,1,time)>= threshold %Comparison with senssivity and thershold
                
                %List of refelections with 8 columns
                detPeaks = [detPeaks; time+eqDistMTX.blockSize/2, ...     %1 timesteps
                    time*eqDistMTX.foreShift/48+eqDistMTX.blockSize/2, ...      %2 time in ms
                    eqDistMTX.mtxData(i,2,1)*180/pi, ...                        %3 azimut in deg
                    eqDistMTX.mtxData(i,3,1)*180/pi, ...                        %4 elevation in deg
                    eqDistMTX.mtxData(i,1,time), ...                            %5 intensity
                    20*log10(eqDistMTX.mtxData(i,1,time)/absoluteMaximum), ...  %6 dB Value
                    eqDistMTX.mtxData(i,2,1), ...                               %7 azimut in rad
                    eqDistMTX.mtxData(i,3,1)];                                  %8 elevation in rad
            end
        end
    end
end

if size(detPeaks,1) < 2
    error('No reflections found. Please adjust parameters!');
end

detPeaks(1,:) = []; %Remove init row

fprintf('\n')


