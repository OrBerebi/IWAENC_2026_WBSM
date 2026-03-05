function Perm=MCRoomPerm(Nmax)

    % this code forms a permute matrix that orders the coefficients that we use
    % to the order MCRoomSim does, following GenSHIndices.m The following does
    % so by C_{MCRoomSIM convention} = Perm x C_{our convention}

    Perm = zeros((Nmax+1)^2);
    sizeP = size(Perm,1);
    ind = 0;
    for n= 0:Nmax

        Perm((ind+1):(ind+1  +(2*n+1) - 1     ),(ind+1):(ind+1  +(2*n+1) - 1     )) = miniMCRoomPerm(n);
        ind = ind + (2*n +1);
    end

    Perm = inv(Perm); 
end