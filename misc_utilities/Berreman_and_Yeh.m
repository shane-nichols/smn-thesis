
kx = 0.1;
ky = 0.2;
epsilon = materialLib('SYEDS', 450);
epsilon = SO3([1,2,3])*diag(epsilon)*SO3([1,2,3]).';

xi = rand(3)/10i;
zeta = -xi.';
mu = eye(3);
%%% METHOD 1
M = [epsilon, xi; zeta, mu]; % 6x6 constitutive matrix
S11 = -M([5 1],[1 5]);
S12 = [-M(5,2),M(5,4);-M(1,2),M(1,4)];
S21 = [M(4,1),M(4,5);-M(2,1),-M(2,5)];
S22 = [M(4,2),-M(4,4);-M(2,2),M(2,4)];
            S31 = [M(6,1)+ky,M(6,5);M(3,1),M(3,5)+kx];
            S32 = [M(6,2)-kx,-M(6,4);M(3,2),-M(3,4)+ky];
            S13 = -[M(5,3)+kx,M(5,6);M(1,3),M(1,6)+ky];
            S23 = [M(4,3)-ky,M(4,6);-M(2,3),-M(2,6)+kx];
S33 = inv(M([6 3],[3 6]));
delta1 = -[S11-S13*S33*S31, S12-S13*S33*S32; S21-S23*S33*S31, S22-S23*S33*S32]
[Psi1,q1] = eig(delta1); % solve for the eigensystem
EHz = [-S33*S31 , -S33*S32] * Psi1;
Efields = [Psi1([1,3],:);EHz(1,:)];

%%% METHOD 2
syms kz         % a symbolic variable for the unknown
K = [0,-kz,ky;kz,0,-kx;-ky,kx,0]; 
q2 = double(solve( det( epsilon + (K + xi)/mu*(K - zeta) ) == 0));
Psi2 = zeros(4);
Efields2 = zeros(3,4);
for i=1:4
    K = [0,-q2(i),ky;q2(i),0,-kx;-ky,kx,0]; % plug in an eigenvalue
    E = null(epsilon + (K + xi)/mu*(K - zeta)); % [x,y,z] componnets of an eigenmodes E-field 
    Efields2(:,i) = E;
    H = mu\(K - zeta)*E; % [x,y,z] componnets of an eigenmodes H-field 
    % next 4 lines make the vector [Ex,Hy,Ey,-Hx]
    Psi2(:,i) = [E(1);H(2);E(2);-H(1)];
end
delta2 = Psi2*diag(q2)/Psi2 % construct delta to show the result matches METHOD 1