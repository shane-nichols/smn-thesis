function Psi = psiAnisoMap(layer, wavelengths, kx_t, n0)

[epsilon, alpha, mu] = materialLib(layer{1}, wavelengths);
eul = layer{3} * pi/180;
Nwl = length(wavelengths);
Ndir = size(kx_t, 2);
Psi = zeros(4,4,Ndir,Nwl);

R = R_ZXZ(eul(1), eul(2), eul(3));
if alpha == 0
    alpha = zeros(3,3,Nwl);
else
    alpha = multiprod(multiprod(R.', 1i .* alpha), R);
end
epsilon = multiprod(multiprod(R.', epsilon), R);
mu = multiprod(multiprod(R.', mu), R);

for i = 1:Nwl
    m1 = epsilon(:,:,i); % views into arrays (faster than repeated indexing)
    m2 = alpha(:,:,i);
    m3 = mu(:,:,i);
    t1 = n0(i);
    for j = 1:Ndir
        R = Rz(kx_t(2, j));
        m4 = R*m2*R.';
        M = [R*m1*R.', -m4; m4.', R*m3*R.'];
        kx = kx_t(1,j) .* t1;
        S11 = -M([5 1],[1 5]);
        S12 = [-M(5,2),M(5,4);-M(1,2),M(1,4)];
        S21 = [M(4,1),M(4,5);-M(2,1),-M(2,5)];
        S22 = [M(4,2),-M(4,4);-M(2,2),M(2,4)];
        S33 = invert2x2(M([6 3],[3 6]));
        S31 = [M(6,1),M(6,5);M(3,1),M(3,5)+kx];
        S32 = [M(6,2)-kx,-M(6,4);M(3,2),-M(3,4)];
        S13 = -[M(5,3)+kx,M(5,6);M(1,3),M(1,6)];
        S23 = [M(4,3),M(4,6);-M(2,3),-M(2,6)+kx];
        t3 = S33*S31;
        t4 = S33*S32;
        delta = -[S11-S13*t3, S12-S13*t4; S21-S23*t3, S22-S23*t4];
        [tPsi, K] = eig(delta, 'vector');
        % eignsystem has to be sorted to use EV matrix as exit medium
        [~, idx] = sort(real(K), 1, 'descend');
        idx = idx([1,2,4,3]);
        Psi(:,:,j,i) = tPsi(:,idx);
    end
end