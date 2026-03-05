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
% W/G/C Wave Generator Core R13-0306
% 
% Copyright (C)2011-2013 Benjamin Bernschütz  
%                        rockzentrale 'AT' me.com
%                        +49 171 4176069 Germany   
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
%  
% [Pnm, kr] = sofia_wgc(N, r, ac, FS, NFFT, AZ, EL, 
%                             t, c, wavetype, ds, lSegLim, uSegLim, SeqN)
% ------------------------------------------------------------------------
% Pnm      Spatial Fourier Coefficients 
%          Columns : nm coeff
%          Rows    : FFT bins 
% kr       kr-Vector 
%          Can also be a matrix [krm; krs] for rigid sphere configurations:
%          [1,:] => krm referring to the Microphone Radius
%          [2,:] => krs referring to the Sphere Radius (Scatterer)
%  
% ------------------------------------------------------------------------   
% N        Maximum transform order
% r        Microphone Radius 
%          Can also be a vector for rigid/dual sphere configurations:
%          [1,1] => rm  Microphone Radius 
%          [2,1] => rs  Sphere Radius or Microphone2 Radius
%          ! If only one radius (rm) is given using a Rigid/Dual Sphere  
%            Configuration: rs = rm and only one kr-vector is returned!
% ac       Array Configuration 
%          0  Open Sphere with p Transducers (NO plc!)
%          1  Open Sphere with pGrad Transducers
%          2  Rigid Sphere with p Transducers
%          3  Rigid Sphere with pGrad Transducers (Thx to Nils Peters!)
%          4  Dual Open Sphere with p Transducers (Thx to Nils Peters!)
% FS       Sampling Frequency
% NFFT     FFT Order (Number of bins) should be 2^x, x=1,2,3,... 
% AZ       Azimuth   angle in [DEG] 0-2pi            
% EL       Elevation angle in [DEG] 0-pi                   
% t        Time Delay in s. The delay has: (t*FS) Samples
% c        Speed of sound in [m/s] (Default: 343m/s)
% wavetype Type of the Wave. 0: Plane Wave (default) 1: Spherical Wave 
% ds       Distance of the source in [m] (For wavetype = 1 only)
%          Warning: If NFFT is smaller than the time the wavefront 
%          needs to travel from the source to the array, the impulse 
%          response will by cyclically shifted (cyclic convolution). 
% ---
% lSegLim  (Lower Segment Limit) Used by the S/W/G wrapper
% uSegLim  (Upper Segment Limit) Used by the S/W/G wrapper
% SegN     (Sement Order)        Used by the S/W/G wrapper
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
