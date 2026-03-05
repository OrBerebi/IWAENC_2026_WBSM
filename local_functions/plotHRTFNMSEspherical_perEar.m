function plotHRTFNMSEspherical_perEar(x, y, omega, freq, freqBands,display_angles,fig_name)
% x, y: [nDirs x nFreqs x 2] estimated and target HRTFs
% omega: [nDirs x 2] = [el, az] in radians
% freq: [nFreqs x 1] frequency vector (Hz)
% freqBands: [nBands x 2] = [f_low, f_high] frequency band ranges

nDirs = size(x,1);
nBands = size(freqBands,1);
nmse_map = zeros(nDirs, nBands, 2);  % [dirs, bands, ear]
mag_map = zeros(nDirs, nBands, 2);  % [dirs, bands, ear]

% Compute NMSE for each direction, frequency band, and ear
for d = 1:nDirs
    for b = 1:nBands
        f_idx = freq >= freqBands(b,1) & freq <= freqBands(b,2);
        for ear = 1:2
            x_band = squeeze(x(d,f_idx,ear));
            y_band = squeeze(y(d,f_idx,ear));
            num = mean(abs(y_band - x_band).^2);
            denom = mean(abs(y_band).^2) + 1e-12;
            nmse_map(d,b,ear) = 10 * log10(num / denom); % in dB
        end
    end
end

% Compute Magnitude error for each direction, frequency band, and ear
for d = 1:nDirs
    for b = 1:nBands
        f_idx = freq >= freqBands(b,1) & freq <= freqBands(b,2);
        for ear = 1:2
            x_band = squeeze(x(d,f_idx,ear));
            y_band = squeeze(y(d,f_idx,ear));
            num = mean(abs(abs(y_band) - abs(x_band)).^2);
            denom = mean(abs(y_band).^2) + 1e-12;
            mag_map(d,b,ear) = 10 * log10(num / denom); % in dB
        end
    end
end

% Convert omega to degrees
az_deg = rad2deg(omega(:,2));                 % [0, 360]
az_deg = mod(az_deg + 180, 360) - 180;
az_deg = -az_deg;

el_deg = rad2deg(pi/2 - omega(:,1));          % [-90, 90]

% Define 2D grid
az_grid = linspace(-179, 180, 360);
el_grid = linspace(-90, 90, 181);
[AZ, EL] = meshgrid(az_grid, el_grid);

figure;
for ear = 1:2
    for b = 1:nBands
        subplot(nBands, 2, 2*(b-1) + ear);

        % Interpolate scattered data to regular grid
        Z = griddata(az_deg, el_deg, nmse_map(:,b,ear), AZ, EL, 'natural');
        %F = scatteredInterpolant(az_deg, el_deg, nmse_map(:,b,ear), 'natural', 'none');
        %Z = F(az_grid, el_grid);
        % Plot as image
        imagesc(az_grid, el_grid, Z);
        set(gca, 'YDir', 'normal');  % North pole on top
        xticks([-90 0 90]);
        yticks([-45 0 45]);
        colormap turbo;
        %xlim([0 display_angles(1)]);
        %ylim([-display_angles(2) display_angles(2)]);
        %ylim([-display_angles(2) display_angles(2)]);
        clim([-15 15]);
        colorbar;
        axis square;
        %title(sprintf('%s Ear: %d–%d Hz', earLabel(ear), freqBands(b,1), freqBands(b,2)));
        title(sprintf('%s Ear: %d-%d Hz', earLabel(ear), freqBands(b,1), freqBands(b,2)), 'Interpreter', 'none');

    end
    xlabel('Azimuth');
    ylabel('Elevation');
end
tmp = fig_name+  "NMSE_sector.fig";
savefig(tmp)

figure;
for ear = 1:2
    for b = 1:nBands
        subplot(nBands, 2, 2*(b-1) + ear);

        % Interpolate scattered data to regular grid
        Z = griddata(az_deg, el_deg, mag_map(:,b,ear), AZ, EL, 'natural');
        %F = scatteredInterpolant(az_deg, el_deg, nmse_map(:,b,ear), 'natural', 'none');
        %Z = F(az_grid, el_grid);
        % Plot as image
        imagesc(az_grid, el_grid, Z);
        set(gca, 'YDir', 'normal');  % North pole on top
        xticks([-90 0 90]);
        yticks([-45 0 45]);
        colormap turbo;
        %xlim([0 display_angles(1)]);
        %ylim([-display_angles(2) display_angles(2)]);
        %ylim([-display_angles(2) display_angles(2)]);
        clim([-25 15]);
        colorbar;
        axis square;
        %title(sprintf('%s Ear: %d–%d Hz', earLabel(ear), freqBands(b,1), freqBands(b,2)));
        title(sprintf('%s Ear: %d-%d Hz', earLabel(ear), freqBands(b,1), freqBands(b,2)), 'Interpreter', 'none');

    end
    xlabel('Azimuth');
    ylabel('Elevation');
end
tmp = fig_name+  "Magnitude_error_sector.fig";
savefig(tmp)

% sgtitle('HRTF NMSE per Direction and Ear (Spherical 2D Projection)');
end

function label = earLabel(ear)
    if ear == 1
        label = 'Left';
    else
        label = 'Right';
    end
end
