function [Psi, n] = psiAmbientInvMap(layer, kx_t, wavelengths)
% aoi : free-space angle of incidence in radians
% layer : layer containing the material name

n = materialLibIso(layer{1}, wavelengths);
a = sqrt(1 - kx_t(1,:).^2); % cos(phi)
a = 1./(2.*a);

Psi = zeros(4, 4, length(a), length(wavelengths));
for j=1:length(wavelengths)
    for k = 1:length(a)
        t1 = 1/(2*n(j));
        t2 = t1 .* a(k) * 2;
        Psi(:,:,k,j) = [a(k),  t1,   0,       0     ;...
                        0 ,      0,   1/2,      t2 ; ...
                        -a(k), t1,   0,        0    ;...
                        0,       0,   1/2,      -t2   ];
    end
end
end