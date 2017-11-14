function layerMatrix = layerBerreman(layer,wavelengths,kx)
% calculates the INVERSE Berreman layer matrix, i.e., the quantity L(-d) in
% J. Opt. Soc. Am. A, 32, 2015, 2049-2057.

[epsilon,alpha,mu] = materialLib(layer{1}, wavelengths);
d = layer{2};
k0 = 2i * pi * d ./ wavelengths;
eul = layer{3} * pi/180;
Nwl = length(wavelengths);
layerMatrix = zeros(4,4,Nwl);

R = R_ZXZ(eul(1), eul(2), eul(3));
if alpha == 0
    alpha = zeros(3,3,Nwl);
else
    alpha = multiprod(multiprod(R.', 1i .* alpha), R);
end

epsilon = multiprod(multiprod(R.', epsilon), R);
mu = multiprod(multiprod(R.', mu), R);

for n = 1:Nwl
    M = [epsilon(:,:,n), -alpha(:,:,n); alpha(:,:,n).', mu(:,:,n)];
    S11 = -M([5 1],[1 5]);
    S12 = [-M(5,2),M(5,4);-M(1,2),M(1,4)];
    S21 = [M(4,1),M(4,5);-M(2,1),-M(2,5)];
    S22 = [M(4,2),-M(4,4);-M(2,2),M(2,4)];
    S33 = [-M(3,6), M(6,6); M(3,3), -M(6,3)] ./ (M(3,3).*M(6,6) - M(3,6).*M(6,3));
    S31 = [M(6,1),M(6,5);M(3,1),M(3,5)+kx(n)];
    S32 = [M(6,2)-kx(n),-M(6,4);M(3,2),-M(3,4)];
    S13 = -[M(5,3)+kx(n),M(5,6);M(1,3),M(1,6)];
    S23 = [M(4,3),M(4,6);-M(2,3),-M(2,6)+kx(n)];
    t1 = S33*S31;
    t2 = S33*S32;
    delta = -[S11-S13*t1, S12-S13*t2; S21-S23*t1, S22-S23*t2];
    [Psi,K] = eig(delta, 'vector');
    layerMatrix(:,:,n) = Psi * diag(exp(K * k0(n))) / Psi;
end
end