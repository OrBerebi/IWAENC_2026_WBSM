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
% Make Matrix 3D Data (HD) R13-0306
%
% Preparing pressure data for visual3d() or mtxToGixel()
% 
% Copyright (C)2011-2013 Benjamin Bernschütz   
%                        rockzentrale 'AT' me.com
%                        +49 171 4176069 Germany  
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
%  
% mtxData = sofia_makeMtx(N, Pnm, dn, krIndex) 
% ------------------------------------------------------------------------     
% mtxData   SOFiA 3D-matrix-data in 1° steps
% ------------------------------------------------------------------------              
% N         Order of the spatial fourier transform     [default = 3]
% Pnm       Spatial Fourier Coefficients (from S/T/C)
% dn        Modal Radial Filters (from M/F)
% krindex   Index of kr Vector                         [default = 1]
% oversize  Integer Factor to increase the resolution. Set oversize = 1
%           (default) to use the mtxData matrix for visual3D(), map3D().
%
% Dependencies: SOFiA P/D/C
% 
% The file generates a SOFiA mtxData Matrix of 181x360 pixels for the
% visualisation with sofia_visual3d() in 1° Steps (65160 plane waves).
% The HD version generally allows to raise the resolution (oversize > 1).  
% (visual3D(), map3D() admit 1° data only, oversize = 1)
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


function mtxData = sofia_makeMTX(N, Pnm, dn, krIndex, oversize)  

disp('SOFiA Visual Matrix Generator Wrapper HD R13-0306');

if nargin < 5
   oversize = 1;
end

oversize = round(oversize);

if oversize < 1 
    oversize = 1;
end
   
angle = zeros(65160*oversize^2,2);
cnt=1;

Elev = linspace(0,180, 181*oversize);
Azim = linspace(0,359, 360*oversize);

for ElevCNT=1:size(Elev,2)
     for AzimCNT=1:size(Azim,2)         
         angle(cnt,:) = [Azim(AzimCNT) Elev(ElevCNT)];
         cnt = cnt+1;
     end
end

angle=angle*pi/180;
Y = sofia_pdc(N, angle, Pnm(:,krIndex), dn(:,krIndex));

cnt=1;
mtxData=zeros(181*oversize,360*oversize);

for ElevCNT=1:size(Elev,2)
    for AzimCNT=1:size(Azim,2)
        mtxData(ElevCNT, AzimCNT)=Y(cnt);
        cnt=cnt+1;
    end
end

