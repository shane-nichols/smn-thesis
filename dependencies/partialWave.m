function MM = partialWave(Psi0, Psi2, layer, wavelengths, kx, bReflect)

%   epsilon and alpha: 3x3 constitutive tensors in the standard setting. 
%   eul: array of ZXZ passive euler angles to rotate the tensors, in radians.
%   AOI: incident angle, in radians

[epsilon,alpha,mu] = materialLib(layer{1}, wavelengths);
d = layer{2};
k0 = 2i * pi * d ./ wavelengths;
eul = layer{3} * pi/180;
Nwl = length(wavelengths);
Ainv = [0.5,0.5,0,0;0,0,0.5,-0.5*1i;0,0,0.5,0.5*1i;0.5,-0.5,0,0];
A = [1,0,0,1;1,0,0,-1;0,1,1,0;0,1i,-1i,0];

R = R_ZXZ(eul(1), eul(2), eul(3));
if alpha == 0
    alpha = zeros(3,3,Nwl);
else
    alpha = multiprod(multiprod(R.', 1i .* alpha), R);
end

epsilon = multiprod(multiprod(R.', epsilon), R);
mu = multiprod(multiprod(R.', mu), R);

Psi1 = zeros(4, 4, Nwl);
Psi1inv = zeros(4, 4, Nwl);
K = zeros(4, Nwl);
for n = 1:Nwl
    M = [epsilon(:,:,n), -alpha(:,:,n); alpha(:,:,n).', mu(:,:,n)];
    S11 = -M([5 1],[1 5]);
    S12 = [-M(5,2),M(5,4);-M(1,2),M(1,4)];
    S21 = [M(4,1),M(4,5);-M(2,1),-M(2,5)];
    S22 = [M(4,2),-M(4,4);-M(2,2),M(2,4)];
    S33 = invert2x2(M([6 3],[3 6]));
    S31 = [M(6,1),M(6,5);M(3,1),M(3,5)+kx(n)];
    S32 = [M(6,2)-kx(n),-M(6,4);M(3,2),-M(3,4)];
    S13 = -[M(5,3)+kx(n),M(5,6);M(1,3),M(1,6)];
    S23 = [M(4,3),M(4,6);-M(2,3),-M(2,6)+kx(n)];
    t1 = S33*S31;
    t2 = S33*S32;
    delta = -[S11-S13*t1, S12-S13*t2; S21-S23*t1, S22-S23*t2];
    [tPsi, tK] = eig(delta, 'vector');
    [~, idx] = sort(real(tK), 1, 'descend');  % sort all eigensystems at once
    idx = idx([1,2,4,3]);
    K(:,n) = tK(idx);
    tPsi1 = tPsi(:,idx);
    Psi1(:,:,n) = tPsi1;
    Psi1inv(:,:,n) = inv(tPsi1);
end
% vectorized Parial Wave calculation
P1 = setDiag(exp(-K(1:2,:) .* k0));
P2 = setDiag(exp( K(3:4,:) .* k0));
Int01 = multiprod(multinv(Psi0), Psi1);
Int12 = multiprod(Psi1inv, Psi2);
Int10 = multiprod(Psi1inv, Psi0);
T01 = invert2x2(Int01([1,2], [1,2], :));
T12 = invert2x2(Int12([1,2],[1,2], :));
R12 = multiprod(Int12([3,4],[1,2], :), T12);
T10 = invert2x2(Int10([3,4],[3,4], :));
R10 = multiprod(Int10([1,2],[3,4], :), T10);
P1 = bigKron(P1);
P2 = bigKron(P2);
T12 = bigKron(T12);
R12 = bigKron(R12);
T10 = bigKron(T10);
R10 = bigKron(R10);

if bReflect
    R01 = multiprod(Int01([3,4], [1,2], :), T01);
    T01 = bigKron(T01);
    R01 = bigKron(R01);
    P2 = multiprod(multiprod(P2,R12),P1); % redefine P2 for in-place 
    MM = R01 + multiprod(multiprod(multiprod(T10, P2), ...
    multinv(eye(4) - multiprod(R10, P2))), T01);
    MM = real(multiprod(multiprod(A, MM), Ainv));
else
    T01 = bigKron(T01);
    MM = multiprod(multiprod(multiprod(T12,P1),...
    multinv(eye(4) - multiprod(multiprod(multiprod(R10,P2),R12),P1))), T01);
    MM = real(multiprod(multiprod(A, MM), Ainv));    
end

end