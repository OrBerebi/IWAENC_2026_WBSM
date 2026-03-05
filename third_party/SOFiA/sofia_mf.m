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
% M/F Modal radial filters R13-0306
%     Soft amplification limiting
%     On-axis powerloss compensation with
%     N0plc to N0 interpolation.
% 
% 
% Copyright (C)2011-2013 Benjamin Bernschütz  
%                        rockzentrale 'at' me.com
%                        +49 171 4176069 Germany  
% 
% Contributions to this version: Nils Peters (nils@icsi.berkeley.edu)
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
% 
% [dn, beam] = SOFIA_MF(N, kr, ac, [a_max], [plc], [fadeover])
% ------------------------------------------------------------------------   
% dn          Vector of modal 0-N frequency domain filters
% beam        Expected free field On-Axis kr-response 
% ------------------------------------------------------------------------
% N           Maximum Order
% kr          Vector or Matrix of kr values
%             First Row   (M=1) N: kr values Microphone Radius
%             Second Row  (M=2) N: kr values Sphere/Microphone2 Radius 
%             [kr_mic;kr_sphere] for Rigid/Dual Sphere Configurations
%             ! If only one kr-vector is given using a Rigid/Dual Sphere  
%             Configuration: kr_sphere = kr_mic 
% ac          Array Configuration: 
%             0  Open Sphere with pressure Transducers (NO plc!)
%             1  Open Sphere with cardioid Transducers
%             2  Rigid Sphere with pressure Transducers
%             3  Rigid Sphere with cardioid Transducers (Thx to Nils Peters!)
%             4  Dual Open Sphere with pressure Transducers (Thx to Nils Peters!)
% a_max       Maximum modal amplification limit in [dB]
% plc         OnAxis powerloss-compensation: 
%             0  Off
%             1  Full kr-spectrum plc
%             2  Low kr only -> set fadeover             
% fadeover    Number of kr values to fade over +/- around min-distance 
%             gap of powerloss compensated filter and normal N0 filters.
%             0 = auto fadeover
% 
 
%
% % CONTACT AND LICENSE INFORMATION:
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
%     Nils Peters                nils 'at' icsi.berkley.edu
% 
