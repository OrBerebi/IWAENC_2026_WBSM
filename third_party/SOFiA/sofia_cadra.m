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
% CAD/R/A CAD Reflection Analysis R13-01-24
%
% Copyright (C)2013 by Maximilian Rühl and Benjamin Bernschütz
%
% This file is part of the SOFiA toolbox under GNU General Public License
%
%
% [] = sofia_cadra(detPeaks, cadPath, offset)
% ------------------------------------------------------------------------
% detPeaks          From sofia_pd()
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
%
% cadPath           Path of the CAD Venue file used by realCADviewer
%
% offset            Free space around the origin. (Distance Compression)
%                   offset [0...1] referring to the maximum distance. 
%                   [Default 0.1] 
%
% This function uses a peak list generated with sofia_pd() and a simple CAD
% model of the respective room to visualize reflections within the CAD
% model.
%
% Dependencies: realCADviewer()



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


function [] = sofia_cadra(detPeaks, cadFile, offset)

if nargin <3, offset        =  0.1; end

if nargin <2
    [File,Path] = uigetfile('','Choose CAD data');
    cadPath=[Path,File];
end

if nargin <1
    error('detPeaks needed');
end

disp('SOFiA CAD/R/A CAD Reflection Analysis R13-01-24');

sphereSizeFactor   = 0.1;
bubblesize         = 50;
sizeOfFont         = 10;

load(cadFile,'venue')

if ~exist('venue','var')
    error('Invalid CAD data')    
end

maxDist = max(max(venue.obj3d(:,3,:)));
maxTime = max(detPeaks(:,1));

[Xm,Ym,Zm] = sph2cart(detPeaks(:,7), pi/2-detPeaks(:,8), detPeaks(:,1)*maxDist*(1-offset)/maxTime+maxDist*offset);

delete(gcf)

[x,y,z] = sphere;
surf(x*sphereSizeFactor, y*sphereSizeFactor, z*sphereSizeFactor) 


set(gcf,'Color', 'w');
axis equal
lighting phong;
camzoom(1.5);
col = colorbar();
xLabel = get(col,'xlabel');
set(xLabel,'String','dB'); 

[Xmb,Ymb,Zmb]=sph2cart(detPeaks(:,7), pi/2-detPeaks(:,8), maxDist);
colormapcopy = colormap;
colorspace = linspace(min(detPeaks(:,6)), max(detPeaks(:,6)), size(colormapcopy,1));
hold on

for i=1:size(detPeaks,1)    
    findMin = abs(colorspace - detPeaks(i,6));
    matchedColor = find(findMin == min(findMin)) ;    
    colorValue = [colormapcopy(matchedColor,1) colormapcopy(matchedColor,2) colormapcopy(matchedColor,3)];    
    line([0 Xmb(i)], [0 Ymb(i)], [0 Zmb(i)],'lineStyle','--','color',colorValue)
    scatter3(Xm(i),Ym(i),Zm(i), bubblesize, detPeaks(i,6), 'marker', 'o', 'MarkerEdgeColor', [0.4 0.4 0.4], 'MarkerFaceColor', colorValue)   
end

hold off

for i=1:size(detPeaks,1)    
    text(Xmb(i),Ymb(i),Zmb(i)+0.2,num2str(i), 'fontsize', sizeOfFont)
end

%3D-Model
realCADviewer(cadFile,'k','k')


