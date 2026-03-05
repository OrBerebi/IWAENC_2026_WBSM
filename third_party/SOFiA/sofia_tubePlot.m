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
% Tube Plot Wrapper R13-0306
%
% Copyright (C)2011-2013 Benjamin Bernschütz
%                        rockzentrale 'AT' me.com
%                        +49 171 4176069 Germany
%
% This file is part of the SOFiA toolbox under GNU General Public License
%
%
% void = sofia_tubePlot(Pnm, dn, FS, cLim, Nmax, frequencyTicks)
% ------------------------------------------------------------------------
% void              No return data.
% ------------------------------------------------------------------------
% Pnm               Spatial Fourier Coefficients (from S/T/C)
% dn                Modal Radial Filters (from M/F)
% cLim              Color scale limits in dB, e.g.: [-40 0] 
%                   Default: Autoscale 
% Nmax              Maximum order (default: max. available in Pnm)
% frequencyTicks    Frequency tick vector in the LOG plot, e.g.: 
%                   [1000, 2000, 4000, 8000], Default: [(100) 1000 10000]
%
% Dependencies: SOFiA P/D/C
%
% The file generates a log plot that depicts the magnitude of the array
% response for a 360° turn in the horizontal plane over the frequency.
% (Origin of the name: Full spectrum plane waves appear as "tubes" in
% the plot.) This plot is extremely helpful for observing the array
% response over the frequency, even if only a single azimutal slice
% is depicted.
%
% Even scenarios that do not match directly can be easily adapted to be
% observable using sofia_tubePlot(). The Fourier coefficients Pnm can
% e.g. be rotated using Wigner-D rotation -> sofia_wdr() in order to be
% in the scope of observation: If you can't turn the telescope, simply
% turn the universe ;)

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


function sofia_tubePlot(Pnm, dn, FS, cLim, Nmax, frequencyTicks)

%  load examplePnm.mat Pnm dn
%  FS = 48000;

Ncoeff = sqrt(size(Pnm,1))-1;
NFFTh  = size(Pnm,2);

if nargin < 5
    Nmax = Ncoeff;
end

Nmax = round(Nmax);

if Nmax > Ncoeff
    Nmax = Ncoeff;
    warning(['The requested order is not delivered in the Pnm coefficients. Reduced to the maximum available order: ',num2str(Nmax)])
end

if Nmax < 0
    Nmax = Ncoeff;
    warning(['Nonsense Order. Set to order: ',num2str(Nmax)])
end

disp('SOFiA Tube Plot Wrapper R13-0306');

radialAngles = [linspace(0,2*pi-(2*pi/360), 360)', pi/2*ones(360,1)];
radialAngles = circshift(radialAngles, 180);
Y            = sofia_pdc(Nmax, radialAngles , Pnm, dn);
scaler       = round(logspace(1,log10(NFFTh),NFFTh));
frequencies  = scaler*(FS/2)/max(scaler);

if nargin < 6
    if size(Pnm,2) > 4000
        frequencyTicks = [100 1000 10000];
    else
        frequencyTicks = [1000 10000];
    end
end

frequencyPointer = 1;
index = 1;

xScaler = [];
while(1)
    if  frequencies(index) >= frequencyTicks(frequencyPointer)
        xScaler(frequencyPointer) = index;
        frequencyPointer = frequencyPointer+1;
    end
    index = index+1;
    if frequencyPointer > size(frequencyTicks,2)
        break
    end
    if index >= size(frequencies,2)
        break
    end
end

logY = zeros(size(Y));
for i=1:size(scaler,2)
    newindex   = scaler(i);
    if newindex == 0
        newindex =1;
    end
    logY(:,i) = Y(:,newindex);
end

set(gcf,'Color','w');

if nargin < 4
    imagesc(20*log10(abs(logY)));
else
    imagesc(20*log10(abs(logY)), cLim);
end

set(gca, 'YDir','normal', 'XTickLabel',[])
set(gca, 'XTICK', xScaler)
set(gca, 'YTICK', [1 90 180 270 359])
set(gca, 'YTickLabel', [char('180'); char('270'); char('  0'); char(' 90'); char('180')])

x_labels = [];

for i=1:length(frequencyTicks)
    x_labels = char(x_labels, num2str(frequencyTicks(i)));
end
x_labels(1,:) = [];

set(gca,'XTickLabel',  x_labels)
xlabel('Frequency in Hz')
ylabel('Azimuth Angle in DEG')

colormap gray

FontSize = 16;
set(gcf, 'Units', 'pixels', 'Position', [100 100 800 380]);
set(gcf, 'PaperPositionMode', 'auto');
set(gca,'LineWidth', 1.5)
set(gcf, 'Color', [1 1 1]); 
set(gca, 'Color', [1 1 1]); 
set(gca, 'FontSize', FontSize);
set(get(gca, 'xlabel'),'FontSize',FontSize)
set(get(gca, 'ylabel'),'FontSize',FontSize)

h = colorbar;
set(get(h,'xlabel'),'String', 'dB');
xlabh = get(h,'xlabel');
set(get(h,'xlabel'),'FontSize', FontSize);
set(xlabh,'Position',get(xlabh,'Position') - [0 .1 0])
set(h,'LineWidth', 1.5);
set(h,'FontSize', FontSize);


