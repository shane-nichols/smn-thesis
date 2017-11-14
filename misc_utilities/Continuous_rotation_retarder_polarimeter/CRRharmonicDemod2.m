function M = CRRharmonicDemod2(f0,f1,p0,p1,d0,d1,I,t)
% demodulation using a harmonic basis of 16 real valued sines/cosines.
% This method is not in my report -- it has similar overhead to LS method,
% but the retardance can be applied after in the linear combinations. 
% The phases must be known to use this minimal basis. Use expanded complex
% basis for applying phases after.

B = zeros(16,1);
A = B;
for index=1:length(t)-1   % demodulation
arg0 = 2*(2*pi*f0.*t(index) + p0);
arg1 = 2*(2*pi*f1.*t(index) + p1);
B(1) = 1/2;
B(2) = sin(arg0);
B(3) = sin(arg1);
B(4) = cos(2*arg0);
B(5) = sin(2*arg0);
B(6) = cos(2*arg1);
B(7) = sin(2*arg1);
B(8) = cos(arg0-arg1);
B(9) = cos(arg0-2*arg1);
B(10) = sin(arg0-2*arg1);
B(11) = cos(2*arg0-arg1);
B(12) = sin(2*arg0-arg1);
B(13) = cos(2*(arg0+arg1));
B(14) = sin(2*(arg0+arg1));
B(15) = cos(2*(arg0-arg1));
B(16) = sin(2*(arg0-arg1));
A = I(index).*B + A;
end
% precalculate stuff for CPU efficiency
M = zeros(4);
a = sin(d0);
b = sin(d1);
c = (1+cos(d0))/2;
d = -(1+cos(d1))/2;
e = (1-cos(d0))/2;
f = (cos(d1)-1)/2;
af = 2/(a*f);
ef = 1/(e*f);
be = 2/(b*e);
sumr = ef*(A(13) + A(15));
sumi = ef*(A(14) + A(16));
difi = ef*(A(14) - A(16));
% construct Mueller matrix
M(1,1) = A(1) - (c/e)*A(4) - (d/f)*A(6) + (c*d)*sumr;
M(1,2) = A(4)/e - d*sumr;
M(1,3) = A(5)/e - d*sumi;
M(2,4) = af*A(10);
M(3,4) = af*A(9);
M(1,4) = A(2)/a - d*M(2,4);
M(2,1) = A(6)/f - c*sumr;
M(3,1) = A(7)/f - c*difi;
M(2,2) = sumr;
M(2,3) = sumi;
M(3,2) = difi;
M(3,3) = -ef*(A(13) - A(15));
M(4,2) = -be*A(12);
M(4,1) = A(3)/b - c*M(4,2);
M(4,3) = be*A(11);
M(4,4) = 2*A(8)/(a*b);
M = M./(length(t)-1)*2;
end

