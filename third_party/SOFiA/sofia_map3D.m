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
% Map3D - Sound Field Visualisation based on Map Projections R12-0104
%
% Copyright (C)2011-2013 Benjamin Bernschütz
%
% This file is part of the SOFiA toolbox under GNU General Public License
%
%
% void = sofia_map3d(mtxData, projection, textSize)
% ------------------------------------------------------------------------
% void         No return.
% ------------------------------------------------------------------------
% mtxData      SOFiA 3D-matrix-data [1°steps]
% projection   Map Projection Type: 
%               'mollweid': Mollweide Projection (#default) 
%               'hammer'  : Hammer Projection
%               'aitoff'  : Hammer/Aitoff Projection
%
% Important: MATLAB Mapping Toolbox Required 
%            http://www.mathworks.de/products/mapping/


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

function sofia_map3D(mtxData, projection)

disp('SOFiA Map3D Visual Map Projection R13-0306');

    mappingToolboxExists = false;
    toolboxes = ver;
    for i = 1: size(toolboxes,2)
         if ~isempty(strfind(toolboxes(i).Name,'Mapping Toolbox'))
             disp('>>> Mapping toolbox: ok');
             mappingToolboxExists = true;
             break
         end
    end
if ~mappingToolboxExists
    disp('>>> Mapping toolbox: failed');
    error('The MATLAB Mapping Toolbox is required and could not be found.');
end

if nargin == 0
   load sofia_exampleMtxData mtxData
end

if nargin < 2
    projection = 'mollweid';
end

if ~strcmp(projection,'mollweid') && ~strcmp(projection,'hammer') && ~strcmp(projection,'aitoff')
error('Type of map projection is not valid');
end

if (size(mtxData,1)~=181 || size(mtxData,2)~=360)
    error('Data type missmatch: SOFiA 3D-matrix-data with 181x360 pixels expected.');
end

mtxData = abs(mtxData);

clf();
set(gcf,'Color', 'w');
mtxData = flipud(mtxData);
mtxData = fliplr(mtxData);
reverse = 1;   % Vertical Label Location Left(0)/Right(1)
textSize = 12; % Size of the Map Labels


% Create referencing matrix
mtxDataR = makerefmat('RasterSize', size(mtxData), ...
    'Latlim', [-90 90], 'Lonlim', [0 360]);

% Set up proj
mapaxes = axesm('MapProjection',projection,'Frame','off','Grid','on','ParallelLabel','on', 'MeridianLabel','on');
if reverse == 0
    setm(gca, 'PLabelMeridian', -90, 'PlineLocation', 30)
else
    setm(gca, 'PLabelMeridian', 90, 'PlineLocation', 30)
end

setm(gca, 'MLabelParallel', 0, 'MlineLocation', 30)

setm(gca, 'FontColor', [0.7, 0.7, 0.7], 'FontSize', textSize, 'FontWeight', 'bold')
setm(gca, 'LabelFormat','none', 'LabelRotation','off', 'maplonlimit', [-180 180]);
spacing = [200 400];
setm(mapaxes, 'glinewidth', 1, 'glinestyle', ':', 'gcolor', [0.7, 0.7, 0.7])

% Display data mapped to the graticule
meshm(mtxData, mtxDataR, spacing);
h=findobj(gcf, 'type', 'text');

if strcmp(projection,'hammer') || strcmp(projection,'aitoff') 
    
    if textSize <= 30
    %Meridian
    set(h(1),'String','');
    set(h(2),'String','');
    set(h(3),'String','30^{\circ}');
    set(h(4),'String','');
    set(h(5),'String','60^{\circ}');
    set(h(6),'String','');
    set(h(7),'String','');
    set(h(8),'String','');
    set(h(9),'String','120^{\circ}');
    set(h(10),'String','');
    set(h(11),'String','150^{\circ}');
    set(h(12),'String','');
    set(h(13),'String','');
    
    %Parallel
    set(h(14),'String','210^{\circ}');
    set(h(15),'String','240^{\circ}');
    set(h(16),'String','270^{\circ}');
    set(h(17),'String','300^{\circ}');
    set(h(18),'String','330^{\circ}');
    
    set(h(20),'String','30^{\circ}');
    set(h(21),'String','60^{\circ}');
    set(h(22),'String','90^{\circ}');
    set(h(23),'String','120^{\circ}');
    set(h(24),'String','150^{\circ}');
    
    else
    %Meridian
    set(h(1),'String','');
    set(h(2),'String','');
   % set(h(3),'String','30^{\circ}');
    set(h(3),'String','');
    set(h(4),'String','');
    set(h(5),'String','');
    set(h(6),'String','');
    set(h(7),'String','');
    set(h(8),'String','');
    set(h(9),'String','');
    set(h(10),'String','');
   % set(h(11),'String','150^{\circ}');
    set(h(11),'String','');
    set(h(12),'String','');
    set(h(13),'String','');
    
    if reverse == 0
        set(h(16),'String','270^{\circ}');
        set(h(22),'String','90^{\circ}');
    else
        %set(h(22),'String','270^{\circ}');
        %set(h(16),'String','90^{\circ}');
        set(h(22),'String','');
        set(h(16),'String','');
        set(h(19),'String','');    
    end
    
    %Parallel
    set(h(14),'String','');
    set(h(15),'String','');    
    set(h(17),'String','');
    set(h(18),'String','');    
    set(h(20),'String','');
    set(h(21),'String','');   
    set(h(23),'String','');
    set(h(24),'String','');
    end
        
elseif strcmp(projection,'mollweid') 
   
    if textSize <= 30   
    %Meridian
    set(h(1),'String','');
    set(h(2),'String','');
    set(h(3),'String','30^{\circ}');
    set(h(4),'String','');
    set(h(5),'String','60^{\circ}');
    set(h(6),'String','');
    set(h(7),'String','');
    set(h(8),'String','');
    set(h(9),'String','120^{\circ}');
    set(h(10),'String','');
    set(h(11),'String','150^{\circ}');
    set(h(12),'String','');
    set(h(13),'String','');
    
    %Parallel
    set(h(14),'String','');
    set(h(15),'String','210^{\circ}');
    set(h(16),'String','240^{\circ}');
    set(h(17),'String','270^{\circ}');
    set(h(18),'String','300^{\circ}');
    set(h(19),'String','330^{\circ}');
    
    set(h(20),'String','0^{\circ}');
    set(h(21),'String','30^{\circ}');
    set(h(22),'String','60^{\circ}');
    set(h(23),'String','90^{\circ}');
    set(h(24),'String','120^{\circ}');
    set(h(25),'String','150^{\circ}');
    set(h(26),'String','');
    
    
    else
    %Meridian
    set(h(1),'String','');
    set(h(2),'String','');
    set(h(3),'String','30^{\circ}');
    set(h(4),'String','');
    set(h(5),'String','');
    set(h(6),'String','');
    set(h(7),'String','');
    set(h(8),'String','');
    set(h(9),'String','');
    set(h(10),'String','');
    set(h(11),'String','150^{\circ}');
    set(h(12),'String','');
    set(h(13),'String','');
    
    %Parallel
    set(h(14),'String','');
    set(h(15),'String','');
    set(h(16),'String','');
    if reverse == 0
        set(h(17),'String','270^{\circ}');
        set(h(23),'String','90^{\circ}');
    else
        set(h(23),'String','270^{\circ}');
        set(h(17),'String','90^{\circ}');
    end
    set(h(20),'String','0^{\circ}');
    set(h(18),'String','');
    set(h(19),'String','');   
    set(h(21),'String','');
    set(h(22),'String','');    
    set(h(24),'String','');
    set(h(25),'String','');
    set(h(26),'String','');        
    end
        
end
axis off
box off
colorbar();



