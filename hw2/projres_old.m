function res = projres(coefs, T, fspace, knodes, zgrid, transmat)

    nz = length(zgrid);
    nk = length(knodes);
    coefs = reshape(coefs,nk,nz);

    res = zeros(nk*nz, 1);

    for iz =1:nz
        z = exp(zgrid(iz));
        znext = exp(zgrid)';

        % C(k,z), l(k,z)
        C = T*coefs(:,iz);
        l = (1./C*.67*z.*knodes.^.33).^(1/1.33);
        % l = ones(nk,1);

        % kp(k,z)
        kp = z*knodes.^.33.*l.^.67 + .9*knodes - C;

        Cnext = zeros(nk,nz);
        for jz = 1:nz
            Cnext(:,jz) = funeval(coefs(:,jz),fspace,kp);
        end
        
        lnext = (1./Cnext*.67.*znext.*kp.^.33).^(1/1.33);
        % lnext = ones(nk,3);
        MRnext = .97 * C./Cnext .* (.33*znext.*kp.^(-.67).*lnext.^.67 + .9);
        res((iz-1)*nk+1:iz*nk) = MRnext * transmat(iz, :)' - 1;
    end

end


