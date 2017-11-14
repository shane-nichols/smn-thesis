function w = BlackmanNuttallPts(n, N)
% returns the nth point of a Blackman-Nuttall window function of N points
a0=0.3635819;
a1=0.4891775;
a2=0.1365995;
a3=0.0106411;
arg = pi*n/N;
w = a0 - a1*cos(2*arg) + a2*cos(4*arg) - a3*cos(6*arg);
end