function out1 = R_ZXZ(t1,t2,t3)
%R_ZXZ
%    OUT1 = R_ZXZ(T1,T2,T3)
%    Active ZXZ Euler rotation matrix

t5 = cos(t3);
t6 = sin(t1);
t7 = cos(t1);
t8 = cos(t2);
t9 = sin(t3);
t10 = sin(t2);
out1 = reshape([t5.*t7-t6.*t8.*t9,-t7.*t9-t5.*t6.*t8,t6.*t10,t5.*t6+t7.*t8.*t9,...
    -t6.*t9+t5.*t7.*t8,-t7.*t10,t9.*t10,t5.*t10,t8],[3,3]);
end
