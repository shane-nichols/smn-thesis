function R = Rz(radians)

c = cos(radians);
s = sin(radians);

R = [c,-s,0;s,c,0;0,0,1];