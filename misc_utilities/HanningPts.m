function w = HanningPts(n,N)
% returns the nth point of a Hanning window function of N points
w = (1-cos(n*2*pi/N))./2;
end
