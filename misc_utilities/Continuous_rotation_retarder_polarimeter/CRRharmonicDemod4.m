function M = CRRharmonicDemod4(f0,f1,p0,p1,d0,d1,I,t)
% Equations 22.
% vectorized (as in MatLab syntax vectorization) complex harmonic demodulation
M = zeros(4);
a = sin(d0);
b = sin(d1);
c = (1+cos(d0))/2;
d = -(1+cos(d1))/2;
e = (1-cos(d0))/2;
f = (cos(d1)-1)/2;

alpha = 2*(2*pi*f0.*t + p0);
beta = 2*(2*pi*f1.*t + p1);
B_1 = sum(I)/2;
B_2 = sum(I.*exp(1i*alpha));
B_3 = sum(I.*exp(1i*beta));
B_4 = sum(I.*exp(2*1i*alpha));
B_5 = sum(I.*exp(2*1i*beta));
B_6 = sum(I.*exp(1i*(alpha-beta)));
B_7 = sum(I.*exp(1i*(alpha-2*beta)));
B_8 = sum(I.*exp(1i*(2*alpha-beta)));
B_9 = sum(I.*exp(1i*2*(alpha+beta)));
B_10 = sum(I.*exp(1i*2*(alpha-beta)));

M(1,1) = real(B_1 - c/e*B_4 - d/f*B_5 + c*d/(e*f)*(B_9 + B_10));
temp = B_4/e - d/(e*f)*(B_9 + B_10);
M(1,2) = real(temp);
M(1,3) = imag(temp);
M(1,4) = imag(B_2/a - 2*d*B_7/(a*f)); 
temp = B_5/f - c/(e*f)*(B_9 + conj(B_10));
M(2,1) = real(temp);
M(3,1) = imag(temp);
temp = (B_9 + B_10)/(e*f);
M(2,2) = real(temp);
M(2,3) = imag(temp);
temp = 2*B_7/(a*f);
M(2,4) = imag(temp);
M(3,4) = real(temp);
temp = (-B_9 + B_10)/(e*f);
M(3,2) = -imag(temp);
M(3,3) = real(temp);
M(4,1) = imag(B_3/b + 2*c*B_8/(b*e));
temp =  2*B_8/(b*e);
M(4,2) = -imag(temp);
M(4,3) = real(temp);
M(4,4) = real(2*B_6/(a*b));
M = M./length(t)*2;
end