function [Psi, n] = psiAmbient(layer, aoi, wavelengths)
% aoi : free-space angle of incidence in radians
% layer : layer containing the material name

n = materialLibIso(layer{1}, wavelengths);
aoi = cos(aoi);

Psi = zeros(4, 4, length(wavelengths), length(aoi));
for k = 1:length(aoi)
    for j = 1:length(wavelengths)
        Psi(:,:,j,k) = [aoi(k),0,-aoi(k),0;n(j),0,n(j),0;0,1,0,1;0,n(j)*aoi(k),0,-n(j)*aoi(k)];
    end
end
end
