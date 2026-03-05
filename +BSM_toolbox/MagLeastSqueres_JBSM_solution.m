function c = MagLeastSqueres_JBSM_solution(V, h, c_priv, W)
    phi = angle(c_priv'*V).';
    H_hat = abs(h).*exp(1j*phi);
    c = (V * V') \ (V * W * conj(H_hat));
end