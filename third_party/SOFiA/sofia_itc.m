% /// ASAR-MARA Research Group
%
% Cologne University of Applied Sciences
% Berlin University of Technology
% University of Rostock
% Deutsche Telekom Laboratories
% WDR Westdeutscher Rundfunk
% IOSONO GmbH
% 
% I/T/C Fast Inverse spatial Fourier Transform Core R13-0306
% 
% Copyright (C)2011-2013 Benjamin Bernschütz  
%                        rockzentrale 'AT' me.com
%                        +49 171 4176069 Germany  
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
% 
% p = sofia_itc(Pnm, angles, [N])
% ------------------------------------------------------------------------
% p      sound pressures (complex data)  
%        Columns : FFT bins
%        Rows    : angles   
% ------------------------------------------------------------------------
% Pnm    spatial Fourier coefficients (e.g. from SOFiA S/T/C)
%        Columns : FFT bins 
%        Rows    : nm coeff  
% 
% angles target angles [AZ1 EL1; AZ2 EL2; ... AZn ELn]
%        Columns : Angle Number 1...n
%        Rows    : AZ EL    
%          
% [N]     *** Optional: Maximum transform order 
%            If not specified the highest order available included in
%            the Pnm coefficients will be taken.
% 
% This is a pure ISFT core that does not involve extrapolation. 
% (=The pressures are referred to the original radius)
%
 
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
