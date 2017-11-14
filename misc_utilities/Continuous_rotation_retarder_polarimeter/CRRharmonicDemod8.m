function M = CRRharmonicDemod8(f0,f1,p0,p1,d0,d1,I,t,polAng1,polAng2)
% Equations 22.
% vectorized (as in MatLab syntax vectorization) complex harmonic demodulation
% This new version generalizes to arbitrary polarizer angles, and also
% demostrates applying the phases after the harmonic demodulation. 
M = zeros(4);
a = -sin(d0);
b = sin(d1);
c = (1+cos(d0))/2;
d = (1+cos(d1))/2;
e = (1-cos(d0))/2;
f = (1-cos(d1))/2;
% demodulation: determine 10 complex harmonic coefficients without retarder
% phases. This is basically a discrete Fourier Transform
alpha = 2*(2*pi*f0.*t);
beta = 2*(2*pi*f1.*t);
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

% apply retarder phases to the coefficients
B_2 = B_2*exp(2i*p0);
B_3 = B_3*exp(2i*p1);
B_4 = B_4*exp(4i*p0);
B_5 = B_5*exp(4i*p1);
B_6 = B_6*exp(2i*(p0-p1));
B_7 = B_7*exp(2i*(p0-2*p1));
B_8 = B_8*exp(2i*(2*p0-p1));
B_9 = B_9*exp(4i*(p0+p1));
B_10 = B_10*exp(4i*(p0-p1));

% precalculate some exponentials of polarizer angles
ep1 = exp(-2i*polAng1);
ep2 = exp(-2i*polAng2);

% take linear combinations of polarizers
M(1,1) = real(B_1 - c/e*ep1^2*B_4 - d/f*ep2^2*B_5 + ...
    c*d/(e*f)*(ep1^2*ep2^2*B_9 + ep1^2*ep2^2'*B_10));
temp = ep1*B_4/e - d/(e*f)*(ep1*ep2^2*B_9 + ep1*ep2^2'*B_10);
M(1,2) = real(temp);
M(1,3) = imag(temp);
M(1,4) = imag(ep1*B_2/a - 2*d*ep1*ep2^2'*B_7/(a*f)); 
temp = ep2*B_5/f - c/(e*f)*(ep1^2*ep2*B_9 + (ep1^2*ep2'*B_10)');
M(2,1) = real(temp);
M(3,1) = imag(temp);
temp = (ep1*ep2*B_9 + ep1*ep2'*B_10)/(e*f);
M(2,2) = real(temp);
M(2,3) = imag(temp);
temp = 2*ep1*ep2'*B_7/(a*f);
M(2,4) = imag(temp);
M(3,4) = real(temp);
temp = (-ep1*ep2*B_9 + ep1*ep2'*B_10)/(e*f);
M(3,2) = -imag(temp);
M(3,3) = real(temp);
M(4,1) = imag(ep2*B_3/b + 2*c*ep1^2*ep2'*B_8/(b*e));
temp =  2*ep1*ep2'*B_8/(b*e);
M(4,2) = -imag(temp);
M(4,3) = real(temp);
M(4,4) = real(2*ep1*ep2'*B_6/(a*b));
M = M./length(t)*8;
end