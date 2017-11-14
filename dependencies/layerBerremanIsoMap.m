function layerMatrix = layerBerremanIsoMap(layer, wavelengths, kx_t, n0)
% calculates the INVERSE Berreman layer matrix for isotropic media. 
%i.e., the quantity L(-d) in J. Opt. Soc. Am. A, 32, 2015, 2049-2057.

% added reshape because some materials give column vectors, which messes
% things up
epsilon = reshape(materialLibIso(layer{1}, wavelengths).^2, size(n0));
Nwl = length(wavelengths);
Ndir = size(kx_t, 2);
d = layer{2};
layerMatrix = zeros(4, 4, Ndir, Nwl);
k0 = (2*pi*d) ./ wavelengths;
for j=1:Nwl
    for k=1:Ndir
        q = sqrt(epsilon(j) - (n0(j).*kx_t(1,k)).^2);
        arg = k0(j) .* q ;
        C = cos(arg);
        S = 1i .* sin(arg);
        qS = q .* S;
        Sq = S ./ q;
        layerMatrix(:,:,k,j) = [...
            C,                    qS ./ epsilon(j),     0,    0  ;...
            epsilon(j) .* Sq,           C,              0,    0  ;...
            0,                          0,              C,    Sq ;...
            0,                          0,              qS,    C ];
    end
end

end
