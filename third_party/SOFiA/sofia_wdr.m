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
% SOFiA W/D/R Wigner-D Rotation R13-0306
%
% Copyright (C)2013 by Benjamin Bernschütz
%
% External routine for the wigner-D matrix calculation: wignerd.p is
% taken from the free EasySpin Toolbox.
% -----------------------------------------------------------------------
% Stefan Stoll and Arthur Schweiger, "EasySpin, a comprehensive software
% package for spectral simulation and analysis in EPR," In: J. Magn.
% Reson. 178(1), 42-55 (2006)
% -----------------------------------------------------------------------
%
% This file is part of the SOFiA toolbox under GNU General Public License
%
% ADVICE: The underlying calculation routine (wigned.p) is NOT part
%         of the SOFiA GNU GPL License!
%
% PnmRot = sofia_wdr(Pnm, xAngle, yAngle, zAngle, deg)
% ------------------------------------------------------------------------
% PnmRot   Output: Rotated spatial Fourier coefficients
%
% Pnm      Input: Spatial Fourier coefficients
% xAngle   Rotation angle around the x-Axis 
% yAngle   Rotation angle around the y-Axis 
% zAngle   Rotation angle around the z-Axis 
% deg      true  (1) - Angles in degree
%          false (0) - Angles in rad (#default if not set)
%
% This routine enables the 3D rotation of the spatial Fourier coefficients
% (Pnm) in the spherical wave spectrum domain.
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
%
%
% The external routine for the wigner-D matrix calculation: wignerd.p is
% taken from the free EasySpin Toolbox (http://www.easyspin.org/)
% -----------------------------------------------------------------------
% Stefan Stoll and Arthur Schweiger, "EasySpin, a comprehensive software
% package for spectral simulation and analysis in EPR," In: J. Magn.
% Reson. 178(1), 42-55 (2006)
% -----------------------------------------------------------------------
% WARNING: The external routine is not part of the SOFiA GNU-GPL License.


function PnmRot = sofia_wdr(Pnm, xAngle, yAngle, zAngle, deg)

disp('SOFiA W/D/R Wigner-D Rotation R13-0306');

if nargin < 5
    deg = false;
end

if deg
    xAngle = xAngle*180/pi;
    yAngle = yAngle*180/pi;
    zAngle = zAngle*180/pi;
end

PnmRot = zeros(size(Pnm));

for i=0:sqrt(size(Pnm,1))-1
    wignerD = wignerd(i,[xAngle, yAngle, zAngle],'-');
    for bin = 1:size(Pnm,2)
        PnmRot(i^2+1:(i+1)^2,bin) = (Pnm(i^2+1:(i+1)^2, bin)' * wignerD')';
    end
end

