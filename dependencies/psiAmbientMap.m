function [Psi, n] = psiAmbientMap(layer, kx_t, wavelengths)

n = materialLibIso(layer{1}, wavelengths);
a = sqrt(1 - kx_t(1,:).^2); % cos(phi)

Psi = zeros(4, 4, length(a), length(wavelengths));
for j=1:length(wavelengths)
    for k = 1:length(a)
        Psi(:,:,k,j) = [a(k),0,-a(k),0;n(j),0,n(j),0;0,1,0,1;0,n(j).*a(k),0,-n(j).*a(k)];
    end
end
end