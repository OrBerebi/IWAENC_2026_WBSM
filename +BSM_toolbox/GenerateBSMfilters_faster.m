function [c_BSM_l, c_BSM_r] = GenerateBSMfilters_faster(BSMobj, V_k, hobj_freq_grid,W)
%% GenerateBSMfilters.m
% Generate the BSM filters - either using complex or magnitude LS
% optimization
%% Inputs:
% BSMobj            : (MATLAB object) containing parameters
% V_k               : (n_mic x Q x freq]) steering vectors (ATF)
% hobj_freq_grid    : HRTF object in earo format interpolated to desired frequencies
%% Outputs:
% c_BSM_l           : (n_mic x freq) BSM filters for left ear in freq. domain
% c_BSM_l           : (n_mic x freq) BSM filters for right ear in freq. domain

    %init
    freqs_sig       = BSMobj.freqs_sig;
    fs = BSMobj.fs;
    magLS           = BSMobj.magLS;
    f_cut_magLS     = BSMobj.f_cut_magLS;
    tol_magLS       = BSMobj.tol_magLS;
    max_iter_magLS  = BSMobj.max_iter_magLS;    
    magLS_cvx       = BSMobj.magLS_cvx;
    n_mic           = BSMobj.n_mic;                    
    normSV          = BSMobj.normSV;
    SNR_lin         = BSMobj.SNR_lin;
    inv_opt         = BSMobj.inv_opt; 
    Tik             = BSMobj.Tik;
    Jan             = BSMobj.Jan;
    lambda          = BSMobj.beta.^2; 
    gamma           = 1e-3;
    omega           = [BSMobj.th_BSMgrid_vec, BSMobj.ph_BSMgrid_vec, BSMobj.r_BSMgrid_vec];
    sector_ang      = BSMobj.sector_ang;
    %
    c_BSM_l     = zeros(n_mic, length(freqs_sig));
    c_BSM_r     = zeros(n_mic, length(freqs_sig)); 
    
    % initialize phase for magLS
    if magLS
        phase_init_l_magLS = pi / 2;    % according to [1]
        phase_init_r_magLS = pi / 2;    % according to [1]        
    end
    
    %normalize steering vectors    
    if normSV
        V_k_norm = vecnorm(V_k, 2, 1);
        V_k = V_k ./ V_k_norm;        
    end
    %
    for f = 1:length(freqs_sig)
        V_k_curr = V_k(:, :, f);        
        %%================= solve for vector c, (V * c = h)
        h_l = hobj_freq_grid.data(:, f, 1);
        h_r = hobj_freq_grid.data(:, f, 2);

        if sum(sum(isnan(V_k_curr)))
            c_BSM_l(:, f) = zeros(n_mic, 1);
            c_BSM_r(:, f) = zeros(n_mic, 1);
        else
            
            % solution to Tikhonov regularization problem            
            if magLS && freqs_sig(f) >= f_cut_magLS                    
                % ===== MAG LS     
                if ~magLS_cvx
                    % Variable exchange method
                    if isequal(W, eye(size(W)))
                        % [c_BSM_l(:, f), phase_last] = BSM_toolbox.TikhonovReg_MagLS_v2(V_k_curr, h_l, (1 / SNR_lin), phase_init_l_magLS, inv_opt, tol_magLS, max_iter_magLS);                                        
                        % phase_init_l_magLS = phase_last;
                        % [c_BSM_r(:, f), phase_last] = BSM_toolbox.TikhonovReg_MagLS_v2(V_k_curr, h_r, (1 / SNR_lin), phase_init_r_magLS, inv_opt, tol_magLS, max_iter_magLS);                                        
                        % phase_init_r_magLS = phase_last;

                        c_BSM_l(:, f) = BSM_toolbox.MagLeastSqueres_BSM_solution(V_k_curr, h_l,c_BSM_l(:, f-1),W, (1 / SNR_lin));
                        c_BSM_r(:, f) = BSM_toolbox.MagLeastSqueres_BSM_solution(V_k_curr, h_r,c_BSM_r(:, f-1),W, (1 / SNR_lin));

                    else
                        if Tik
                            c_BSM_l(:, f) = BSM_toolbox.MagLS_Tik_BSM_solution(V_k_curr, h_l,lambda,gamma,omega,sector_ang,c_BSM_l(:, f-1));
                            c_BSM_r(:, f) = BSM_toolbox.MagLS_Tik_BSM_solution(V_k_curr, h_r,lambda,gamma,omega,sector_ang,c_BSM_r(:, f-1));
                        elseif Jan
                            c_BSM_l(:, f) = BSM_toolbox.MagLeastSqueres_JBSM_solution(V_k_curr, h_l,c_BSM_l(:, f-1),W);
                            c_BSM_r(:, f) = BSM_toolbox.MagLeastSqueres_JBSM_solution(V_k_curr, h_r,c_BSM_r(:, f-1),W);    
                        else
                            c_BSM_l(:, f) = BSM_toolbox.MagLeastSqueres_BSM_solution(V_k_curr, h_l,c_BSM_l(:, f-1),W, (1 / SNR_lin));
                            c_BSM_r(:, f) = BSM_toolbox.MagLeastSqueres_BSM_solution(V_k_curr, h_r,c_BSM_r(:, f-1),W, (1 / SNR_lin));    
                        end
                    end
                else
                    % SDP with CVX toolbox
                    c_BSM_l(:, f) = BSM_toolbox.TikhonovReg_MagLS_CVX(V_k_curr, h_l, (1 / SNR_lin), inv_opt);
                    c_BSM_r(:, f) = BSM_toolbox.TikhonovReg_MagLS_CVX(V_k_curr, h_r, (1 / SNR_lin), inv_opt);
                    % ==============
                end


            else
                % ===== CMPLX LS 
                if Tik 
                    c_BSM_l(:, f) = BSM_toolbox.LeastSqueres_Tik_BSM_solution(V_k_curr, h_l,lambda,gamma,omega,sector_ang);
                    c_BSM_r(:, f) = BSM_toolbox.LeastSqueres_Tik_BSM_solution(V_k_curr, h_r,lambda,gamma,omega,sector_ang);
                elseif Jan
                    c_BSM_l(:, f) = BSM_toolbox.LeastSqueres_JBSM_solution(V_k_curr, h_l,W);
                    c_BSM_r(:, f) = BSM_toolbox.LeastSqueres_JBSM_solution(V_k_curr, h_r,W);       
                else
                    if isequal(W, eye(size(W)))
                        %c_BSM_l(:, f) = BSM_toolbox.TikhonovReg(V_k_curr, h_l, (1 / SNR_lin), inv_opt);
                        %c_BSM_r(:, f) = BSM_toolbox.TikhonovReg(V_k_curr, h_r, (1 / SNR_lin), inv_opt);

                        c_BSM_l(:, f) = BSM_toolbox.LeastSqueres_BSM_solution(V_k_curr, h_l,W, (1 / SNR_lin));
                        c_BSM_r(:, f) = BSM_toolbox.LeastSqueres_BSM_solution(V_k_curr, h_r,W, (1 / SNR_lin));
                    else
                        c_BSM_l(:, f) = BSM_toolbox.LeastSqueres_BSM_solution(V_k_curr, h_l,W, (1 / SNR_lin));
                        c_BSM_r(:, f) = BSM_toolbox.LeastSqueres_BSM_solution(V_k_curr, h_r,W, (1 / SNR_lin));       
                    end
                    

                    
                end
            end
        end            
    end 

    % c_BSM_f_cmplx = cat(3,c_BSM_l,c_BSM_r);
    % nFFT = 2*(size(c_BSM_f_cmplx,2)-1);
    % c_BSM_f_full = zeros(size(c_BSM_f_cmplx,1), nFFT, size(c_BSM_f_cmplx,3));
    % c_BSM_f_full(:, 1:nFFT/2+1, :) = c_BSM_f_cmplx;
    % c_BSM_f_full(:, nFFT/2+2:end, :) = conj(flip(c_BSM_f_cmplx(:, 2:nFFT/2, :), 2));
    % 
    % n_shift = nFFT/2;
    % k = (0:nFFT-1)';  % frequency bin indices
    % shift_factor = exp(-1j * 2 * pi * k * n_shift / nFFT);  % column vector
    % c_BSM_f_full = permute(c_BSM_f_full,[2,1,3]);
    % c_BSM_f_full = c_BSM_f_full .* shift_factor;  % automatically broadcast across channels
    % c_BSM_f_full = permute(c_BSM_f_full,[2,1,3]);
    % c_BSM_f_full = c_BSM_f_full(:,1:nFFT/2+1,:);
    % c_BSM_l      = c_BSM_f_full(:,:,1);
    % c_BSM_r      = c_BSM_f_full(:,:,2);
    % 
    % 
    % c_BSM_f = cat(3,c_BSM_l,c_BSM_r);
    % nFFT = 2*(size(c_BSM_f,2)-1);
    % c_BSM_f(:, end+1 : nFFT,:) = 0;
    % c_BSM_t = ifft(c_BSM_f, [], 2, 'symmetric');
    % c_BSM_t_zf = zero_phase_it(c_BSM_t,nFFT,fs);
    % c_BSM_f_zf = fft(c_BSM_t_zf, [], 2);
    % c_BSM_f_zf = c_BSM_f_zf(:,1:nFFT/2+1,:);
    % c_BSM_l      = c_BSM_f_zf(:,:,1);
    % c_BSM_r      = c_BSM_f_zf(:,:,2);





    
end




