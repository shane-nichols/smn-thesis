function [Psi, n] = psiAmbientInv(layer, aoi, wavelengths)
% aoi : free-space angle of incidence in radians
% layer : layer containing the material name

n = materialLibIso(layer{1}, wavelengths);
aoi = 1./(2.*cos(aoi));

Psi = zeros(4, 4, length(wavelengths), length(aoi));
for k = 1:length(aoi)
    for j = 1:length(wavelengths)
        t1 = 1/(2*n(j));
        t2 = t1 .* aoi(k) * 2;
        Psi(:,:,j,k) = [aoi(k),  t1,   0,       0     ;...
                        0 ,      0,   1/2,      t2 ; ...
                        -aoi(k), t1,   0,        0    ;...
                        0,       0,   1/2,      -t2   ];
    end
end
end