function plotSourcesOnSphere(Omega, title_name, array_mic_positons)
% plotSourcesOnSphere - Plots source directions on a unit sphere
%
% Syntax:
%   plotSourcesOnSphere(Omega, title_name, array_mic_positons)
%
% Input:
%   Omega - [DOA x 2] matrix of [elevation, azimuth] angles in radians
%           elevation: angle from z-axis (0 = north pole), [0, pi]
%           azimuth: angle in xy-plane from x-axis, [0, 2*pi]
%   array_mic_positons - [M x 2] matrix of mic positions [elevation, azimuth]

    % Convert spherical coordinates (elevation, azimuth) to Cartesian
    elevation = Omega(:, 1);
    azimuth   = Omega(:, 2);
    mic_el = array_mic_positons(:,1);
    mic_az = array_mic_positons(:,2);

    ear_el = [pi/2 ; pi/2];
    ear_az = [pi/2 ; 3*pi/2];

    x = sin(elevation) .* cos(azimuth);
    y = sin(elevation) .* sin(azimuth);
    z = cos(elevation);

    x_mic = sin(mic_el) .* cos(mic_az);
    y_mic = sin(mic_el) .* sin(mic_az);
    z_mic = cos(mic_el);

    x_ear = sin(ear_el) .* cos(ear_az);
    y_ear = sin(ear_el) .* sin(ear_az);
    z_ear = cos(ear_el);
    
    % Define colors using MATLAB's "lines" colormap
    cmap = lines(3); % Get 3 distinguishable colors
    source_color = cmap(1, :);
    mic_color    = cmap(2, :);
    ear_color    = cmap(3, :);
    
    % Plot the unit sphere
    [sx, sy, sz] = sphere(50); % Higher resolution
    figure;
    surf(sx, sy, sz, 'FaceAlpha', 0.1, 'EdgeColor', 'none'); 
    hold on;
    axis equal;
    
    % Plot points with appropriate sizes and colors
    h_sources = scatter3(x, y, z, 50, 'filled', 'MarkerFaceColor', source_color);
    h_mics    = scatter3(x_mic, y_mic, z_mic, 100, 'filled', 'MarkerFaceColor', mic_color);
    h_ears    = scatter3(x_ear, y_ear, z_ear, 200, 'filled', 'MarkerFaceColor', ear_color);
    
    % Axes labels and title
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    %title(title_name);
    grid on;
    view(135, 30);
    
    % Add legend
    legend([h_sources, h_mics, h_ears], {'Sources', 'Microphones', 'Ears'}, 'Location', 'best');
end
