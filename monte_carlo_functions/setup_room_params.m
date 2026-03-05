function roomsetup = setup_room_params(type_str,diffuese_flag,specular_flag)
    switch type_str
        case 'small'
            roomsetup.type       = type_str;
            roomsetup.roomDim    = [6, 4, 2.5];
            roomsetup.T60_target = 0.2;
        case 'medium'
            roomsetup.type       = type_str;
            roomsetup.roomDim    = [10, 7, 3.5];
            roomsetup.T60_target = 0.4;
        case 'large'
            roomsetup.type       = type_str;
            roomsetup.roomDim    = [20, 15, 6];
            roomsetup.T60_target = 0.8;
    end
    
    poss_coeff = 0.31;
    roomsetup.arrayPos = [roomsetup.roomDim(1)*poss_coeff, roomsetup.roomDim(2)*poss_coeff, 1.7];

    % --- FIX: Calculate Reflection Coefficient (Amplitude) ---
    % Standard ISM Simulators expect Reflection Amplitude (r), not Absorption (alpha).
    % Eyring Relation: r = exp( - (0.161 * V) / (2 * S * T60) )
    
    L = roomsetup.roomDim(1); W = roomsetup.roomDim(2); H = roomsetup.roomDim(3);
    V = L*W*H;
    S = 2*(L*W + L*H + W*H);
    target_T60 = roomsetup.T60_target;

    % Note the factor of 2 in denominator. 
    % Without it, you get Energy Reflection (R^2). With it, you get Amplitude (R).
    R_amplitude_theoretical = exp( - (0.161 * V) / (2 * S * target_T60) ); 
    
    CorrectionFactor = 1; 
    roomsetup.R = R_amplitude_theoretical * CorrectionFactor;
    
    % --- FIX: Pass Reflection Amplitude Directly ---
    % critical_distance_diffuse computes (1 - R^2), so it expects R to be Amplitude.
    % Previously, passing (1-R) was incorrect because R was Absorption.
    roomsetup.critical_distance = RoomParams.critical_distance_diffuse(roomsetup.roomDim, roomsetup.R);

    % Rest of setup
    roomsetup.S          = 0.1; 
    roomsetup.freq_bands = [125, 250, 500, 1000, 2000, 4000];
    
    roomsetup.N_PW       = 18;
    roomsetup.fs         = 16e3;
    %roomsetup.fs         = 48e3;
    %roomsetup.N_PW       = 32;
    roomsetup.roomsim_fs = 48e3;
    roomsetup.c          = 343;

    [roomsetup.max_dist, roomsetup.limit_azimuth] = get_safe_dist_and_angle(roomsetup);

    roomsetup.shutup          = true;
    roomsetup.direct_flag     = true;
    roomsetup.diffuese_flag   = diffuese_flag;
    roomsetup.specular_flag   = specular_flag;
    roomsetup.sig_gender      = 'omnidirectional';
end

function [max_dist, limit_azimuth] = get_safe_dist_and_angle(roomsetup)
    x = roomsetup.arrayPos(1);
    y = roomsetup.arrayPos(2);
    L = roomsetup.roomDim(1);
    W = roomsetup.roomDim(2);
    
    d_right = L - x;
    d_front = W - y;
    d_left = x;
    d_back = y;

    all_dists  = [d_right, d_front, d_left, d_back];
    all_angles = [0,       90,      180,    -90];

    [max_dist, idx] = min(all_dists);
    limit_azimuth = all_angles(idx);
end