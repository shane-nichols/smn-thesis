function CRRlsTEST(f0,f1,p0,p1,d0,d1,t)
arg0 = 2*(2*pi*f0.*t + p0);
arg1 = 2*(2*pi*f1.*t + p1);
Cp0 = cos(arg0);
C2p0 = cos(2*arg0);
Sp0 = sin(arg0);
S2p0 = sin(2*arg0);
Cp1 = cos(arg1);
C2p1 = cos(2*arg1);
Sp1 = sin(arg1);
S2p1 = sin(2*arg1);
Cd0 = cos(d0);
Sd0 = sin(d0);
Cd1 = cos(d1);
Sd1 = sin(d1);
a = sin(d0);
b = sin(d1);
c = (1+cos(d0))/2;
d = -(1+cos(d1))/2;
e = (1-cos(d0))/2;
f = (cos(d1)-1)/2;

%M(1,1) = (0.5-(d*e*C2p1+c*C2p0.*(f-2*d*C2p1))/(e*f))*I;
D_A(1) = sqrt(2)*(f-2*d*C2p1)/(2*f);
D_A(2) = sqrt(2)*C2p1/f;
D_A(3) = sqrt(2)*S2p1/f;
D_A(4) = sqrt(2)*Sp1/f;
D_G(1) = sqrt(2)*(e-2*c*C2p0)/(2*e);
D_G(2) = sqrt(2)*C2p0/e;
D_G(3) = sqrt(2)*S2p0/e;
D_G(4) = sqrt(2)*Sp0/e;
Dmatrix = inv(transpose(D_A)*D_A)
inv(Dmatrix);
B_G(1) = 1;
B_G(2) = Cp0.^2 + Cd0.*Sp0.^2;
B_G(3) = Cp0.*Sp0.*(1 - Cd0);
B_G(4) = Sd0.*Sp0;
B_A(1) = 1;
B_A(2) = -(Cp1.^2 + Cd1.*Sp1.^2);
B_A(3) = -Cp1.*Sp1.*(1 - Cd1);
B_A(4) = Sd1.*Sp1;
(1/(kron(B_A,B_G)*kron(B_A,B_G)'))*kron(B_A,B_G);



end


