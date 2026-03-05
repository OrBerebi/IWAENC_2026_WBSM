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
% S/F/E Sound Field Extrapolation R13-0306
% 
% Copyright (C)2011-2013 Benjamin Bernschütz 
%                        rockzentrale 'AT' me.com
%                        Germany +49 171 4176069  
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
% Pnm_krb = sofia_sfe(Pnm_kra, kra, krb, problem) 
% ------------------------------------------------------------------------     
% Pnm_krb Spatial Fourier Coefficients, extrapolated to rb
% ------------------------------------------------------------------------              
% 
% Pnm_kra Spatial Fourier Coefficients from SOFiA S/T/C
% 
% kra     k*ra Vector
% krb     k*rb Vector
% problem 1: Interior problem [default]  2: Exterior problem

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

function Pnm_krb = sofia_sfe(Pnm_kra, kra, krb, problem) 

disp('SOFiA S/F/E - Sound Field Extrapolation R13-0306');
disp(' ');

if nargin < 3
    error('Arguments missing. 3 Args required: Pnm_krb = sofia_sfe(Pnm_kra, kra, krb)')
end

if nargin < 4
   problem = 1;
end

if size(kra,2) ~= size(krb,2) ||  size(kra,2) ~= size(Pnm_kra,2) 
   error('FFT bin number or dimension mismatch. Pnm_kra, kra and krb must have the same M-dimension.') 
end

FCoeff  = size(Pnm_kra,1);
N       = sqrt(FCoeff)-1;

nvector=zeros(FCoeff,1);

index = 1;
for n=0:N
    for m=-n:n
        nvector(index) = n;
        index = index+1;
    end
end

nvector = repmat(nvector,1,size(Pnm_kra,2));
kra     = repmat(kra,FCoeff,1);
krb     = repmat(krb,FCoeff,1);

if problem == 2
    hn_kra  = sqrt(pi./(2*kra)).*besselh(nvector+.5,1,kra);
    hn_krb  = sqrt(pi./(2*krb)).*besselh(nvector+.5,1,krb);
    exp     = hn_krb./hn_kra;
else
    jn_kra  = sqrt(pi./(2*kra)).*besselj(n+.5,kra);
    jn_krb  = sqrt(pi./(2*krb)).*besselj(n+.5,krb);
    plot(abs(jn_kra'))
    exp     = jn_krb./jn_kra;
    if ~isempty(find(abs(exp)>1e2)) %40dB
       disp('WARNING: Extrapolation might be unstable for one or more frequencies/orders!');
    end
end

Pnm_krb = Pnm_kra.*exp;

