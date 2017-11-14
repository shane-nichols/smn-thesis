function layerMatrix = layerBerremanIso(layer, wavelengths, kx)
% calculates the INVERSE Berreman layer matrix for isotropic media. 
%i.e., the quantity L(-d) in J. Opt. Soc. Am. A, 32, 2015, 2049-2057.

% added reshape because some materials give column vectors, which messes
% things up
epsilon = reshape(materialLibIso(layer{1}, wavelengths).^2, size(kx));

d = layer{2};
layerMatrix = zeros(4,4,length(wavelengths));

q = sqrt(epsilon-kx.^2);
arg = (2*pi*d) .* (q ./ wavelengths);
C = cos(arg);
S = 1i .* sin(arg);
qS = q .* S;
Sq = S ./ q;

layerMatrix(1,1,:) = C;
layerMatrix(1,2,:) = qS ./ epsilon;
layerMatrix(2,1,:) = epsilon .* Sq;
layerMatrix(2,2,:) = C;
layerMatrix(3,3,:) = C;
layerMatrix(3,4,:) = Sq;
layerMatrix(4,3,:) = qS;
layerMatrix(4,4,:) = C;

end
