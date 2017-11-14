function R = Rx(radians)

c = cos(radians);
s = sin(radians);

R = [1,0,0;0,c,-s;0,s,c];