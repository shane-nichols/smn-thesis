function M = CRRlsDemodNew(f0,f1,p0,p1,d0,d1,I,t)
% Equation 18.
% Analytic least-Squares demodulation.  
% It can be optimized more than this, it's written to be readable.

arg0 = 2*(2*pi*f0.*t + p0);
arg1 = 2*(2*pi*f1.*t + p1);

C2p0 = cos(2*arg0);
Sp0 = sin(arg0);
S2p0 = sin(2*arg0);
C2p1 = cos(2*arg1);
Sp1 = sin(arg1);
S2p1 = sin(2*arg1);

a = -sin(d0);
b = -sin(d1);
c = (1+cos(d0))/2;
d = -(1+cos(d1))/2;
e = (1-cos(d0))/2;
f = (cos(d1)-1)/2;

DG1 = 1/2 - C2p0*c/e;
DG2 = C2p0/e;
DG3 = S2p0/e;
DG4 = Sp0/a;
DA1 = 1/2 - C2p1*d/f;
DA2 = C2p1/f;
DA3 = S2p1/f;
DA4 = Sp1/b;

M(1,1) = (DA1.*DG1)*I;
M(1,2) = (DA1.*DG2)*I;
M(1,3) = (DA1.*DG3)*I;
M(1,4) = (DA1.*DG4)*I;
M(2,1) = (DA2.*DG1)*I;
M(2,2) = (DA2.*DG2)*I;
M(2,3) = (DA2.*DG3)*I;
M(2,4) = (DA2.*DG4)*I;
M(3,1) = (DA3.*DG1)*I;
M(3,2) = (DA3.*DG2)*I;
M(3,3) = (DA3.*DG3)*I;
M(3,4) = (DA3.*DG4)*I;
M(4,1) = (DA4.*DG1)*I;
M(4,2) = (DA4.*DG2)*I;
M(4,3) = (DA4.*DG3)*I;
M(4,4) = (DA4.*DG4)*I;
M = 16*M./length(t);

end