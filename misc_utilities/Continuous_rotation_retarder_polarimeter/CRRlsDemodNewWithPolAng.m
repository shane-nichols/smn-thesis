function M = CRRlsDemodNewWithPolAng(f0,f1,p0,p1,d0,d1,I,t,polAng1,polAng2)
% Equation 18.
% Analytic least-Squares demodulation.  
% It can be optimized more than this, it's written to be readable.

arg0 = 2*(2*pi*f0.*t + p0);
arg1 = 2*(2*pi*f1.*t + p1);

a = -sin(d0);
b = sin(d1);
c = (1+cos(d0))/2;
d = (1+cos(d1))/2;
e = (1-cos(d0))/2;
f = (1-cos(d1))/2;

DG1 = 1/2 - cos(2*arg0-4*polAng1)*c/e;
DG2 = cos(2*arg0-2*polAng1)/e;
DG3 = sin(2*arg0-2*polAng1)/e;
DG4 = sin(arg0-2*polAng1)/a;
DA1 = 1/2 - cos(2*arg1-4*polAng2)*d/f;
DA2 = cos(2*arg1-2*polAng2)/f;
DA3 = sin(2*arg1-2*polAng2)/f;
DA4 = sin(arg1-2*polAng2)/b;

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