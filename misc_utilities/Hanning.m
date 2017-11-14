function w = Hanning(N)
% returns a Hanning window of N points
n = 0:1:(N-1);
arg = 2*pi/N;
w = (1-cos(n*arg))./2;
end
