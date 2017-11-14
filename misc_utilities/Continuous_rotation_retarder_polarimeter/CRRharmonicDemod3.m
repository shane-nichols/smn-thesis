function M = CRRharmonicDemod3(f0,f1,p0,p1,d0,d1,I,t)
% vectorized version of CRRharmonicDemod2 - much faster
% for some reason it runs a bit slower than lsDemod. Perhaps it can be
% optimized more because it should be similar. 

% demodulation
arg0 = 2*(2*pi*f0.*t + p0);
arg1 = 2*(2*pi*f1.*t + p1);
A1 = sum(I)/2;
A2 = sum(I.*sin(arg0));
A3 = sum(I.*sin(arg1));
A4 = sum(I.*cos(2*arg0));
A5 = sum(I.*sin(2*arg0));
A6 = sum(I.*cos(2*arg1));
A7 = sum(I.*sin(2*arg1));
A8 = sum(I.*cos(arg0-arg1));
A9 = sum(I.*cos(arg0-2*arg1));
A10 = sum(I.*sin(arg0-2*arg1));
A11 = sum(I.*cos(2*arg0-arg1));
A12 = sum(I.*sin(2*arg0-arg1));
A13 = sum(I.*cos(2*(arg0+arg1)));
A14 = sum(I.*sin(2*(arg0+arg1)));
A15 = sum(I.*cos(2*(arg0-arg1)));
A16 = sum(I.*sin(2*(arg0-arg1)));

% precalculate stuff for efficiency
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
sumr = ef*(A13 + A15);
sumi = ef*(A14 + A16);
difi = ef*(A14 - A16);
% construct Mueller matrix
M(1,1) = A1 - (c/e)*A4 - (d/f)*A6 + (c*d)*sumr;
M(1,2) = A4/e - d*sumr;
M(1,3) = A5/e - d*sumi;
M(2,4) = af*A10;
M(3,4) = af*A9;
M(1,4) = A2/a - d*M(2,4);
M(2,1) = A6/f - c*sumr;
M(3,1) = A7/f - c*difi;
M(2,2) = sumr;
M(2,3) = sumi;
M(3,2) = difi;
M(3,3) = -ef*(A13 - A15);
M(4,2) = -be*A12;
M(4,1) = A3/b - c*M(4,2);
M(4,3) = be*A11;
M(4,4) = 2*A8/(a*b);
M = M./length(t)*8;
end

