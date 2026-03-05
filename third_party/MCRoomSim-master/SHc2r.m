function Perm=SHc2r(Nmax)

    % this code forms a permute matrix from the Normalized Complex Spherical Harmonics to
    % the Normalized Real Spherical Harmonics
    % Perm matrix hold the relation- Ynm_{Real} = Perm x Ynm_{Complex}

    Perm = zeros((Nmax+1)^2);
    sizeP = size(Perm,1);
    ind = 0;
    for n= 0:Nmax

        Perm((ind+1):(ind+1+(2*n+1)-1),(ind+1):(ind+1+(2*n+1)-1)) = miniSHc2r(n);
        ind = ind + (2*n +1);
    end

    Perm=conj(Perm);
   
end