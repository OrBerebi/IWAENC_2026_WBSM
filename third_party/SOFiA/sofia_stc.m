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
% S/T/C Fast Spatial Fourier Transform Core R13-0306
% 
% Copyright (C)2011-2013 Benjamin Bernschütz  
%                        rockzentrale 'AT' me.com
%                        +49 171 4176069 Germany  
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
%  
% Pnm = sofia_stc(N, fftData, grid)
% ------------------------------------------------------------------------
% Pnm      Spatial Fourier Coefficients 
%          Columns : nm coeff
%          Rows    : FFT bins   
% ------------------------------------------------------------------------   
% N        Maximum transform order
% 
% fftdata  Frequency domain sounfield data, e.g. from sofia_fdt()) 
%          Columns : number of spatial sampling position 
%          Rows    : FFT bins (complex sound pressure data)
% 
%          ! WARNING Cast fftdata to double if core crashes:
%            fftData = cast(fftData ,'double');
% 
% grid     Sample grid configuration
%          Columns : s=1...S spatial positions
%          Rows    : [AZ_s EL_s GridWeight_s]
%          AZ in [0...2pi] and EL [0...pi] in RAD
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
