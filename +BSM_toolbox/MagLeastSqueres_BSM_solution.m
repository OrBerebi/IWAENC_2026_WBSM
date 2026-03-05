function c = MagLeastSqueres_BSM_solution(V, h, c_priv, W,lambda)
    phi = angle(c_priv'*V).';
    H_hat = abs(h).*exp(1j*phi);
    % Solve with Tikhonov regularization
    c = (V * abs(W).^2 * V' + lambda * eye(size(V,1))) \ (V * abs(W).^2 * conj(H_hat));
end