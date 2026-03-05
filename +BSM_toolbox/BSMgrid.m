%% BSMgrid.m
% Returns the directions of assumed sources in BSM 
%% Inputs:
% source_distribution : (scalar) 0 - nearly uniform (t-design), 1 - spiral nearly uniform
% Q                   : (scalar) Assumed number of sources
%% Outputs:
% th_BSMgrid_vec      :  (Q x 1) elevation angles 
% ph_BSMgrid_vec      :  (Q x 1) azimuth angles 

function [th_BSMgrid_vec, ph_BSMgrid_vec] = BSMgrid(source_distribution, Q, angles)

    switch source_distribution
        case 0
            % calculate steering vectors with nearly uniform sampling
            N_NU_sources_direction = min( [floor(sqrt(Q)) - 1, 10] );
            %[~, th_BSMgrid_vec, ph_BSMgrid_vec] = BSM_toolbox.uniform_sampling_extended(N_NU_sources_direction);
            [~, th_BSMgrid_vec, ph_BSMgrid_vec] = uniform_sampling_extended(N_NU_sources_direction);
            th_BSMgrid_vec = th_BSMgrid_vec.';
            ph_BSMgrid_vec = ph_BSMgrid_vec.';
            ph_BSMgrid_vec = ph_BSMgrid_vec + pi;
            
            % Value of pi (180 degrees in radians)
            pi_value = pi;
            
            % Tolerance for ±1 degree in radians
            tolerance = deg2rad(0.5);
            
            % Logical index for values within ±1 degree of pi
            index_to_remove = abs(th_BSMgrid_vec - pi_value) <= tolerance;
            
            % Remove elements from both vectors
            th_BSMgrid_vec(index_to_remove) = [];
            ph_BSMgrid_vec(index_to_remove) = [];

            %plot_sampling(th_BSMgrid_vec, ph_BSMgrid_vec);
            %plot_sampling_hammer(th_BSMgrid_vec, ph_BSMgrid_vec);
        case 1
            % calculate steering vectors with spiral nearly uniform sampling
            [~, th_BSMgrid_vec, ph_BSMgrid_vec] = SpiralSampleSphere_w_weights(Q);

            %plot_sampling(th_grid_vec, ph_grid_vec);
            %plot_sampling_hammer(th_grid_vec, ph_grid_vec);
        case 2 %Lebedev Sampaling
            lebedev_order = Q;
            [gridData_Inf, ~, ~] = sofia_lebedev(Order_to_degree_lebedev(lebedev_order),0);
            
            th_BSMgrid_vec = gridData_Inf(:,2);
            ph_BSMgrid_vec = gridData_Inf(:,1);

            %plot_sampling(th_BSMgrid_vec, ph_BSMgrid_vec);
        case 3
            % this is a sector filtred lebedev grid
            lebedev_order = 35;
            [gridData_Inf, ~, ~] = sofia_lebedev(Order_to_degree_lebedev(lebedev_order),0);
            th_BSMgrid_vec = gridData_Inf(:,2);
            ph_BSMgrid_vec = gridData_Inf(:,1);
            r_BSMgrid_vec = ones(size(th_BSMgrid_vec));
            omega = [th_BSMgrid_vec, ph_BSMgrid_vec, r_BSMgrid_vec];
           
            
            [omega_fov,~] = filterAzEl(omega, angles);
            th_BSMgrid_vec = omega_fov(:,1);
            ph_BSMgrid_vec = omega_fov(:,2);


            % th_BSMgrid_vec = omega(:,1);
            % ph_BSMgrid_vec = omega(:,2);
        case 4
            % non filtred lebedev
            lebedev_order = 35;
            [gridData_Inf, ~, ~] = sofia_lebedev(Order_to_degree_lebedev(lebedev_order),0);
            th_BSMgrid_vec = gridData_Inf(:,2);
            ph_BSMgrid_vec = gridData_Inf(:,1);
            r_BSMgrid_vec = ones(size(th_BSMgrid_vec));
            omega = [th_BSMgrid_vec, ph_BSMgrid_vec, r_BSMgrid_vec];
            th_BSMgrid_vec = omega(:,1);
            ph_BSMgrid_vec = omega(:,2);


    end
end





