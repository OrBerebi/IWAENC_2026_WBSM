function c = LeastSqueres_BSM_solution(V, h,W,lambda)
    % ORIGINAL  
    %c = pinv(V*V')*V*conj(h);

    % Wighted
    c = (V * abs(W).^2 * V'+ lambda * eye(size(V,1))) \ (V * abs(W).^2 * conj(h));


    % [U, S, Vh] = svd(V, 'econ');  % V ≈ U*S*Vh'
    % 
    % s = diag(S);
    % tol = max(size(V)) * eps(max(s));
    % r = sum(s > tol);
    % 
    % U_r = U(:, 1:r);
    % S_r = S(1:r, 1:r);
    % 
    % Vin_h = V * conj(h);
    % c = U_r * (S_r \ (S_r \ (U_r' * Vin_h)));  % S_r^{-2} part

    


    % cond_num = cond(G);
    % disp(num2str(cond_num));
    % if cond_num > 1e3
    %     warning('Condition number is high: %.2e. The solution may be unstable.', cond_num);
    % end

    % Compute Gram matrix
    % G = V * V';
    % Compute least-squares solution
    % c = pinv(G) * V * conj(h);

end