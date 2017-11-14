function psi = psiIsoMap(layer, wavelengths, kx_t, n0)

epsilon = materialLibIso(layer{1},wavelengths).^2;
Ndir = size(kx_t, 2);
Nwl = length(wavelengths);
psi = zeros(4, 4, Ndir, Nwl);


% the vectorized code below is equivalent to this
for j = 1:Nwl
    for k = 1:Ndir
        q = sqrt(epsilon(j) - (n0(j) .* kx_t(1,k)).^2);
        t1 = q/epsilon(j);
        t2 = 1/q;
        psi(:,:,j) = [   t1,           0,          0,         -t2;...
                          1,           0,          0,           1;...
                          0,          t2,         -t1,          0;...
                          0,           1,          1,           0];
    end
end


end
    
    