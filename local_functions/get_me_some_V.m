function [V_k, V_t] = get_me_some_V(BSMobj,Omega)

    BSMobj.normSV = true;
    [V_k,V_t] = BSM_toolbox.CalculateSteeringVectors_v2(BSMobj, BSMobj.N_PW,Omega(:,1), Omega(:,2)); % [freq x Omega x mics]


    % tmp = V_k;
    % tmp(end+1 : BSMobj.filt_samp, :, :) = 0;
    % V_t = ifft(tmp, [], 1, 'symmetric');

    %I WANNA CHANGED THIS TO be definced by the array radius, smaple rate
    %and the speed of sound -> to calcualte the larges posible sample delay
    %max_delay = 2* round((BSMobj.r_array(1)./BSMobj.c).*BSMobj.fs);
    %V_t = circshift(V_t, max_delay, 1);
    
    V_k = permute(V_k, [2 1 3]); % [Omega x freq x mics]
    V_t = permute(V_t, [2 1 3]); % [Omega x time x mics]

    V_k = permute(V_k, [3 1 2]); % [mics x Omega x freq]
    V_t = permute(V_t, [3 1 2]); % [mics x Omega x time]

    
end