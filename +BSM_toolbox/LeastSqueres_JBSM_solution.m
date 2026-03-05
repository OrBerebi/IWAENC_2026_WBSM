function c = LeastSqueres_JBSM_solution(V, h,W)
    
    c = (V * V') \ (V * W * conj(h));

end