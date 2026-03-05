function c = LeastSqueres_Tik_BSM_solution(V,h,lambda,gamma,omega,sector_ang)
    [~,mask_FoV] = filterAzEl(omega, sector_ang);
    Vs = V(:,mask_FoV);
    Vo = V(:,~mask_FoV);

    hs = h(mask_FoV);
    ho = h(~mask_FoV);

    M1 = pinv(Vs*Vs' + lambda*(Vo*Vo')+gamma*eye(size(Vo,1)));
    M2 = Vs*conj(hs)+lambda*Vo*conj(ho);
    c = M1*M2;
end