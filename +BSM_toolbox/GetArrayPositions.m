function [th_array, ph_array, ph_rot_array,r_array] = GetArrayPositions(arrayType, n_mic, array_rot_az)
%% GetArrayPositions.m
% Calculate array positions (spherical coordinates) given arrayType
%%Inputs:
% arrayType           : (scalar) 0 - spherical array, 1 - semi-circular array, 2 - full-circular array
% n_mic               : (scalar) number of mics
% array_rot_az        : (scalar) degree of azimuth rotation (rad)
%%Outputs:
% th_array            : (1 x n_mic) Microphone elevation
% ph_array            : (1 x n_mic) Original Microphone azimuth
% ph_rot_array        : (1 x n_mic) Rotated Microphone azimuth

    if arrayType == 0
        % %================= generate coordinates of spherical open array        
        % [~, th_array, ph_array] = SpiralSampleSphere_w_weights(n_mic);
        % th_array = th_array.';
        % ph_array = ph_array.'; 

        %================= generate coordinates of spherical Lebedev array
        %(Order 1 = 6 mics)
        % LebedevOrderAvaliable   = [ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65];
        % Mic2Lebedev_convert = [0,0,0,0,0,1,0,0,0,0,0,0,0,2];
        % degrees_avail=[6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230,... 
        %     266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702,... 
        %     3074, 3470, 3890, 4334, 4802, 5294, 5810];
        % 
        % lebedev_order           = find(degrees_avail==n_mic);
        % [gridData_Inf, ~, ~]    = sofia_lebedev(Order_to_degree_lebedev(lebedev_order),0);
        % 
        % ph_array = gridData_Inf(:,1).'; 
        % th_array = gridData_Inf(:,2).';

        % N_NU_sources_direction = min( [floor(sqrt(n_mic)) - 2, 10] );
        N_uniform_sampling_mic_number = [4,12,32,36,73,84,108,144];
        [val,N_uniform_sampling_order] = min(abs(N_uniform_sampling_mic_number - n_mic));
        [~, th_array, ph_array] = uniform_sampling_extended(N_uniform_sampling_order);

        r_array = 0.1* ones(size(th_array));

    elseif arrayType == 1
        % parameters of semi circular array
        n_mic_tmp = n_mic +2; 

        mic_idx = 1:n_mic_tmp;
        ph_array = deg2rad(linspace(-90,90,n_mic_tmp));
        %ph_array = deg2rad(90 - 180 / (n_mic - 1) * (mic_idx - 1));
        ph_array = wrapTo2Pi(ph_array);
        ph_array = ph_array(2:end-1);
        th_array = deg2rad(90) * ones(size(ph_array));
        r_array = 0.1* ones(size(th_array));

    elseif arrayType == 2
        
        rot_resolution = floor(360/n_mic);
        start_ang = 0;
        ph_array = deg2rad(start_ang:rot_resolution:360-rot_resolution);
        %ph_array = deg2rad(start_ang:rot_resolution:360-start_ang);
        ph_array = wrapTo2Pi(ph_array);
        th_array = deg2rad(90) * ones(size(ph_array));
        r_array = 0.1* ones(size(th_array));
            
    elseif arrayType == 3

        mic_L = [-90,75,-20]*1e-3;


        mic1 = [-29,82,-5]*1e-3;
        mic2 = [30,-1,-1]*1e-3;
        mic3 = [11,-77,-2]*1e-3;
        mic4 = [-60,-83,-5]*1e-3;

        mic_R = [-90,-75,-20]*1e-3;

        mics_pos_car = [mic_L;mic1;mic2;mic3;mic4;mic_R];

        % Fix the mics positions in space (Translation)
        ears_y = [0.09,-0.09];
        ears_x = [0,0];
        ears_z = [0,0];
        v_l = [ears_x(1),ears_y(1),ears_z(1)];
        v_r = [ears_x(2),ears_y(2),ears_z(2)];
        v_1 = [mics_pos_car(1,1),mics_pos_car(1,2),mics_pos_car(1,3)];
        v_6 = [mics_pos_car(6,1),mics_pos_car(6,2),mics_pos_car(6,3)];
        t_l = v_l - v_1;
        t_r = v_r - v_6;
        vec_t = [t_l;t_r];
        vec_t = mean(vec_t,1);

        mics_pos_car(:,1) = mics_pos_car(:,1)+ vec_t(1);
        mics_pos_car(:,2) = mics_pos_car(:,2)+ vec_t(2);
        mics_pos_car(:,3) = mics_pos_car(:,3)+ vec_t(3);


        [th_array,ph_array,r_array] = c2s(mics_pos_car);
    
    elseif arrayType == 4

        ph_array = deg2rad([53.3544,-0.4775,-37.3210,-70.1278 ]);
        th_array = deg2rad([81.6503,81.0032,81.9333,80.3541]);
        r_array = 0.08* ones(size(th_array));
    elseif arrayType == 5
        mic_L = [-90,75,-20]*1e-3;

        mic1 = [-29,82,-5]*1e-3;
        mic2 = [30,-1,-1]*1e-3;
        mic3 = [11,-77,-2]*1e-3;

        mic_R = [-90,-75,-20]*1e-3;

        mics_pos_car = [mic2;mic1;mic3;mic_L;mic_R];

        % Fix the mics positions in space (Translation)
        ears_y = [0.09,-0.09];
        ears_x = [0,0];
        ears_z = [0,0];
        v_l = [ears_x(1),ears_y(1),ears_z(1)];
        v_r = [ears_x(2),ears_y(2),ears_z(2)];
        v_1 = [mics_pos_car(4,1),mics_pos_car(4,2),mics_pos_car(4,3)];
        v_6 = [mics_pos_car(5,1),mics_pos_car(5,2),mics_pos_car(5,3)];
        t_l = v_l - v_1;
        t_r = v_r - v_6;
        vec_t = [t_l;t_r];
        vec_t = mean(vec_t,1);

        mics_pos_car(:,1) = mics_pos_car(:,1)+ vec_t(1);
        mics_pos_car(:,2) = mics_pos_car(:,2)+ vec_t(2);
        mics_pos_car(:,3) = mics_pos_car(:,3)+ vec_t(3);


        [th_array,ph_array,r_array] = c2s(mics_pos_car);


    elseif arrayType == 6
        %mics_pos_car = 1e-2*[9.95,-4.76,0.68;10.59, 0.74, 5.07;9.95, 4.49, 0.76; 9.28, 6.41, 5.12; 9.93, -5.66, 5.22; -0.42, -8.45, 3.35; -0.48, 7.75, 3.49];


        mics_pos_car = 1e-2*[9.95,-4.76,0.68;10.59, 0.74, 5.07;9.95, 4.49, 0.76; 9.28, 6.41, 5.12; 9.93, -5.66, 5.22];

        [th_array,ph_array,r_array] = c2s(mics_pos_car);
    end
    % rotate array in azimuth
    ph_rot_array = wrapTo2Pi(ph_array + 0); 
    %ph_rot_array = wrapTo2Pi(ph_array + array_rot_az);

end



