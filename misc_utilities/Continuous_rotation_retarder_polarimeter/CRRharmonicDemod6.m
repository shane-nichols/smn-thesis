function M = CRRharmonicDemod6(f0,f1,p0,p1,d0,d1,I,t)
% matrix version of harmonic demodulation. not efficient, I only made this so I
% could output the big matrix to LaTex format after testing it. But then, I
% never put this into derivation into my report, so.. 
% 
% demodulation
arg0 = 2*(2*pi*f0.*t + p0);
arg1 = 2*(2*pi*f1.*t + p1);
A1 = sum(I);
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

A = [A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16]';

M = zeros(4);
a = sin(d0);
b = sin(d1);
c = (1+cos(d0))/2;
d = -(1+cos(d1))/2;
e = (1-cos(d0))/2;
f = (cos(d1)-1)/2;

R(1,1)=1/2;
R(1,4)=-c/e;
R(1,6)=-d/f;
R(1,13)=c*d/(e*f);
R(1,15)=c*d/(e*f);
R(2,4)=1/e;
R(2,13)=-d/(e*f);
R(2,15)=-d/(e*f);
R(3,5)=1/e;
R(3,14)=-d/(e*f);
R(3,16)=-d/(e*f);
R(4,2)=1/a;
R(4,10) = -2*d/(a*f);
R(5,6) = 1/f;
R(5,13) = -c/(e*f);
R(5,15) = -c/(e*f);
R(6,13) = 1/(e*f);
R(6,15) = 1/(e*f);
R(7,14) = 1/(e*f);
R(7,16) = 1/(e*f);
R(8,10) = 2/(a*f);
R(9,7) = 1/f;
R(9,14) = -c/(e*f);
R(9,16) = c/(e*f);
R(10,14) = 1/(e*f);
R(10,16) = -1/(e*f);
R(11,13) = -1/(e*f);
R(11,15) = 1/(e*f);
R(12,9) = 2/(a*f);
R(13,3) = 1/b;
R(13,12) = 2*c/(b*e);
R(14,12) = -2/(b*e);
R(15,11) = 2/(b*e);
R(16,8) = 2/(a*b);
M = reshape(R*A,4,4)';
M = M./length(t)*2;
end


