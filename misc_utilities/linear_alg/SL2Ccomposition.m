av = 3*randn(1,3);
bv = 3*randn(1,3);
%a = sqrt(sum(av.^2));
a = sqrt(sum(av.^2));
b = sqrt(sum(bv.^2));
ea = av./a;
eb = bv./b;
Pauli_a = (pauli(2)*ea(1) + pauli(3)*ea(2) + pauli(4)*ea(3));
Pauli_b = (pauli(2)*eb(1) + pauli(3)*eb(2) + pauli(4)*eb(3));
expm(-a*1i*Pauli_a/2);
G1 = cos(a/2)*eye(2) - 1i*sin(a/2)*Pauli_a;
expm(-b*1i*Pauli_b/2);
G2 = cos(b/2)*eye(2) - 1i*sin(b/2)*Pauli_b;

K=1/sqrt(det(G1));
T = 2*acos(K*(G1(1,1)+G1(2,2))/2);
O = K*T/(2*sin(T/2));
p_out = [1i*2*log(1/K),1i*O*(G1(1,1)-G1(2,2)),1i*O*(G1(1,2)+G1(2,1)),O*(G1(2,1)-G1(1,2))].';

G2*G1
abv = (sin(a/2)*cos(b/2).*ea + cos(a/2)*sin(b/2).*eb + ...
    sin(a/2)*sin(b/2).*cross(eb,ea) );
Pauli_ab = (pauli(2)*abv(1) + pauli(3)*abv(2) + pauli(4)*abv(3));

(cos(a/2)*cos(b/2) - sin(a/2)*sin(b/2)*(ea*eb.'))*eye(2) - 1i*Pauli_ab;
(cos(a/2)*cos(b/2) - sin(a/2)*sin(b/2)*(ea*eb.'))*eye(2) - 1i*Pauli_ab;

% ab = sqrt(sum(abv.^2));
% eab = abv./ab;
% Pauli_eab = (pauli(2)*eab(1) + pauli(3)*eab(2) + pauli(4)*eab(3));
% cos(ab/2)*eye(2) - 1i*sin(ab/2)*Pauli_eab
 ab = 2*acos(cos(a/2)*cos(b/2) - sin(a/2)*sin(b/2)*(ea*eb.'))
 (a-b) - ea*eb.';
 ab2 = 2*asin(sqrt(sum(abv.^2)))
 eab = abv./sqrt(sum(abv.^2));
 
Pauli_eab = (pauli(2)*eab(1) + pauli(3)*eab(2) + pauli(4)*eab(3));
cos(ab/2)*eye(2) - 1i*sin(ab/2)*Pauli_eab

sl2c(eab*ab)




