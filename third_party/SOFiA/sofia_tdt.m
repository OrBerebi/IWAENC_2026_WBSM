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
% T/D/T - Time Domain Transform R13-0306
% 
% Copyright (C)2011-2013 Benjamin Bernschütz   
% 
% This file is part of the SOFiA toolbox under GNU General Public License
% 
%  
% timeDomainSignal = sofia_tdt(Y, [win], [resampleFactor], [minPhase])
% ------------------------------------------------------------------------     
% timeDomainSignal   Reconstructed Time Domain Signal
%                    Columns : Index / Channel: IR1, IR2, ... IRn
%                    Rows    : Impulse responses (time domain)
% ------------------------------------------------------------------------              
% Y                  Frequency domain FFT data for multiple channels
%                    Columns : Index / Channel
%                    Rows    : FFT data (frequency domain)
% 
% [win]              Window Signal tail [0...1] with a HANN window
%                    0    off (#default)
%                    0-1  window coverage (1 full, 0 off)
%                   
% [resampleFactor]   Optional resampling: Resampling factor 
%                    e.g. FS_target/FS_source
%                    Resampling is done using MATLAB RESAMPLE 
%                    (See MATLAB documentation for more details)
%                    ! Signal Processing Toolbox required  
%
% [minPhase]         Optional minimum phase reduction 
%                    0 off (#default)
%                    1 on
%                    ! Signal Processing Toolbox required 
%
% 
% This function recombines time domain signals for multiple channels from
% frequency domain data. It is made to work with half-sided spectrum FFT 
% data.  The impulse responses can be windowed.  The IFFT blocklength is 
% determined by the Y data itself:
% 
% Y should have a size [NumberOfChannels x ((2^n)/2)+1] with n=[1,2,3,...] 
% and the function returns [NumberOfChannels x resampleFactor*2^n] samples.
% 
% Dependencies: MATLAB Signal Processing Toolbox required for 
%               windowing and resampling.    
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



function timeDomainSignal = sofia_tdt(Y, win, resampleFactor, minPhase)

disp('SOFiA T/D/T - Time Domain Transform R13-0306');

if nargin == 0
   error('Arguments missing: timeDomainSignal = sofia_tdt(Y, [win], [resampleFactor])');
end

if nargin < 2
   win = 0;
end

if win > 1
   error('Argument win must be in the range 0 to 1.')
end

if nargin < 4
   minPhase = 0;
end

signalProcessingExists = 0;
toolboxes = ver;
for i = 1: size(toolboxes,2)
    if ~isempty(strfind(toolboxes(i).Name,'Signal Processing Toolbox'))
       signalProcessingExists = 1;   
       break
    end
end

if ~signalProcessingExists 
   disp('WARNING: Signal processing toolbox not found. Resampling and minimum phase disabled.') 
end

Y(:,1) = Y(:,2); 
y = [Y(:,:), conj(fliplr(Y(:,2:end-1)))];
y = real(ifft(y,[],2)); 

if signalProcessingExists && minPhase ~= 0    
    y    = [y, zeros(size(y))]';
    Y    = fft(y);
    Y(Y == 0) = 1e-21;
    img  = imag(hilbert(log(abs(Y))));
    y    = real(ifft(abs(Y) .* exp(-1i*img)));
    y    = y(1:end/2,:)';    
end

if win ~= 0   
    irLength = size(y,2);
    j = floor(irLength*win):irLength;
    winfkt = 0.5+0.5*cos(2*pi*(j-((irLength-1)/2))/(irLength-1)); 
    y(:,end-size(winfkt,2)+1:end) = y(:,end-size(winfkt,2)+1:end).*repmat(winfkt,size(y,1),1);
end

if nargin == 3 && signalProcessingExists 
    timeDomainSignal = resample(y',resampleFactor,1)';
else
    timeDomainSignal = y;   
end

