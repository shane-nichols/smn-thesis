function M = CRRlsDemod2(f0,f1,p0,p1,d0,d1,I,t)
% Equation 18.
% Analytic least-Squares demodulation.  
% It can be optimized more than this, it's written to be readable.

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
a = -sin(d0);
b = -sin(d1);
c = (1+cos(d0))/2;
d = -(1+cos(d1))/2;
e = (1-cos(d0))/2;
f = (cos(d1)-1)/2;

%M(1,1) = (0.5-(d*e*C2p1+c*C2p0.*(f-2*d*C2p1))/(e*f))*I;
M(1,1) = ((e-2*c*C2p0).*(f-2*d*C2p1))/(2*e*f)*I;
M(1,2) = C2p0.*(f-2*d*C2p1)/(e*f)*I;
M(1,3) = (f-2*d*C2p1).*S2p0/(e*f)*I;
M(1,4) = (f-2*d*C2p1).*Sp0/(a*f)*I;
M(2,1) = (e-2*c*C2p0).*C2p1/(e*f)*I;
M(2,2) = 2*C2p0.*C2p1/(e*f)*I;
M(2,3) = 2*C2p1.*S2p0/(e*f)*I;
M(2,4) = 2*C2p1.*Sp0/(a*f)*I;
M(3,1) = (e-2*c*C2p0).*S2p1/(e*f)*I;
M(3,2) = 2*C2p0.*S2p1/(e*f)*I;
M(3,3) = 8*Cp0.*Cp1.*Sp0.*Sp1/(e*f)*I;
M(3,4) = 4*Cp1.*Sp0.*Sp1/(a*f)*I;
M(4,1) = (e-2*c*C2p0).*Sp1/(b*e)*I;
M(4,2) = 2*C2p0.*Sp1/(b*e)*I;
M(4,3) = 4*Cp0.*Sp0.*Sp1/(b*e)*I;
M(4,4) = 2/(a*b)*Sp0.*Sp1*I;
M = M./length(t)*2;
end