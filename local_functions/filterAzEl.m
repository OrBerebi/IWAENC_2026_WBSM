function [omega_fov,mask] = filterAzEl(omega, angles)
    % angles = [el_target_deg, az_target_deg]
    el_target = deg2rad(angles(1));
    az_target = deg2rad(angles(2));
    
    % Convert zenith angle (theta) to elevation (el = pi/2 - theta)
    omega_el = pi/2 - omega(:,1);  % Now: omega(:,1) = elevation
    % omega(:,2) = azimuth
    
    % Create mask for azimuth in sector around ±az_target
    mask_az = abs(wrapToPi(omega(:,2))) <= az_target;
    
    % Create mask for elevation in sector around ±el_target
    mask_el = abs(omega_el) <= el_target;
    
    % Combine masks
    mask = mask_az & mask_el;
    
    % Apply mask
    omega_fov     = omega(mask,:);

end