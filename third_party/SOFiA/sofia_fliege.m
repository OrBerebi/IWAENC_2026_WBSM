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
% SOFiA Fliege Grid R14-0131
% 
% Copyright (C)2011-2014 by Benjamin Bernschütz 
%
% This file is part of the SOFiA toolbox under GNU General Public License
%
%  
% [gridData, Npoints, Nmax] = sofia_lebede(degree, plot)
% ------------------------------------------------------------------------     
%
% gridData           Fliege quadrature including weigths(W):
%                    [AZ_1 EL_1 W_1;
%                     AZ_2 EL_2 W_2;
%                     ...
%                     AZ_n EL_n W_n]
%
% Npoints            Total number of nodes
% Nmax               Highest stable grid order  
%
% ------------------------------------------------------------------------
% 
% Order              Fliege Degree (Number of nodes)
%                    Call sofia_fliege() to obtain a 
%                    list of valid degrees.
%
% plot               Show a globe plot of the selected grid 
%                    0: Off, 1: On [default]
% 
% This function computes Lebedev quadrature nodes and weigths
% in the SOFiA/VariSphear data format.
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
% Copyright (C)2011-2014 Benjamin Bernschütz [1,2] et al.(§)
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
% External routines for the quadrature calculation in: lebedev_calc.m
% -----------------------------------------------------------------------
% Author: Rob Parrish, The Sherrill Group, CCMST Georgia Tech
% Email : robparrish 'at' gmail.com
% Date  : 03/24/2010
% -----------------------------------------------------------------------
% WARNING: The external routine is not part of the SOFiA GNU-GPL License.


function [gridData, Npoints, Nmax] = sofia_fliege(degree, plot)

disp('SOFiA Fliege Grid R14-0131');

degrees_avail = [4 9 16 25 36 49 64 81 100 121 144 169 196 225 256 289 324 361 400 441 484 529 576 625 676 729 784 841 900];

if nargin<1
   show_valid_degrees();
   gridData =[];
   Npoints  =[];
   Nmax     =[];
   return
end

if nargin<2
    plot = true;
end

if isempty(find(degrees_avail==degree,1))
   error('Invalid quadrature degree. Call sofia_fliege() to obtain a list of valid degrees.');
end

try
    load sofia_FliegeQuads.mat
catch
    disp('Quadrature data file >> sofia_FliegeQuads.mat << not accessible!');
end

varName = ['Fliege_',num2str(degree)];   
eval(['gridData = ',varName,';']);

Npoints = size(gridData,1);
Nmax    = floor(sqrt(size(gridData,1))-1);

if plot
    plot_grid(gridData);
end
end

function show_valid_degrees()
disp(' ');
disp('Valid degrees (number or nodes) are: ');
disp(' ');
disp('4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225, 256, 289, ')
disp('324, 361,400, 441, 484, 529, 576, 625, 676, 729, 784, 841, 900')
end

function plot_grid(gridData)

[Xm,Ym,Zm]=sph2cart(gridData(:,1),gridData(:,2)-pi/2,1.01);

colormap Gray;

if size(Xm,1)>1500
    plot3(Xm,Ym,Zm,'marker','.','markerfacecolor','g','color','g','linestyle','none')
else
    plot3(Xm,Ym,Zm,'marker','o','markerfacecolor','g','color','g','linestyle','none')
end
axis off;
hold on;
grid off;
sphere;
axis equal;
rotate3d on;
light;
alpha(.8);
lighting phong;
camzoom(1.4);
hold off;
end

