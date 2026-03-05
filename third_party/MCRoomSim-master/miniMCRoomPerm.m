function perm=miniMCRoomPerm(n)

    % a help function for MCRooPerm, permuting for each given n.

    perm = zeros((2*n+1));
    sizeP = size(perm,1);
    perm((floor(sizeP/2)+1),(2*n+1)) = 1;
    for ii= 1:(floor(sizeP/2))
        perm((floor(sizeP/2)+1-ii),(2*n+1) - 2*ii +1 ) = 1;
        perm((floor(sizeP/2)+1+ii),(2*n+1) - 2*ii ) = 1;
    end
end