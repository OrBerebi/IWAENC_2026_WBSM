function c_BSM = get_BSM_filters_CVX(V_k_bsm, hobj_bsm, omega_bsm, omega_bsm_FoV, epsilon,lambda)
%c_l    = SoftConstrainedEpsLS_CVX(V_k_bsm, hobj_bsm.data(:,:,1), omega_bsm, omega_bsm_FoV, epsilon,lambda);
%c_r    = SoftConstrainedEpsLS_CVX(V_k_bsm, hobj_bsm.data(:,:,2), omega_bsm, omega_bsm_FoV, epsilon,lambda);
tic;
c_l    = constrainedEpsLS_CVX(V_k_bsm, hobj_bsm.data(:,:,1), omega_bsm, omega_bsm_FoV, epsilon);
c_r    = constrainedEpsLS_CVX(V_k_bsm, hobj_bsm.data(:,:,2), omega_bsm, omega_bsm_FoV, epsilon);
c_out = cat(3,c_l,c_r);
c_BSM.f_ls    = c_out;

c_l    = constrainedEpsMagLS_CVX_v2(V_k_bsm, hobj_bsm.data(:,:,1), omega_bsm, omega_bsm_FoV, 0);
c_r    = constrainedEpsMagLS_CVX_v2(V_k_bsm, hobj_bsm.data(:,:,2), omega_bsm, omega_bsm_FoV, 0);

c_out = cat(3,c_l,c_r);
c_BSM.f_mls    = c_out;

elapsed_time = toc;
fprintf('CVX execution time: %.6f seconds\n', elapsed_time);

end