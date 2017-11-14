function R = Ry(radians)

c = cos(radians);
s = sin(radians);

R = [c,0,s;0,1,0;-s,0,c];