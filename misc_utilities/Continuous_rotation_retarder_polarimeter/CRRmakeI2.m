function I = CRRmakeI2(f1,f2,p1,p2,d1,d2,M,t,polAng1,polAng2)
% simulation of the continuous rotating retarder light intensity. Here, we
% are able to set the polarizer angles to any value.

% I = array of light intensity values at detector
% f0 = frequency of first retarder in Hz
% f1 = frequency of second retarder in Hz
% p0 = phase of first retarder in radians
% p1 = phase of second retarder in radians
% d0 = retardance of first retarder in radians
% d1 = retardance of second retarder in radians
% M = 4x4 test Mueller matrix
% t = time values to compute I
% polAng0 = angle of the input polarizer, in radians
% polAng1 = angle of the output polarizer, in radians

polAng1 = polAng1*2;
polAng2 = polAng2*2;
Cd1 = cos(d1);
Sd1 = sin(d1);
Cd2 = cos(d2);
Sd2 = sin(d2);
P1 = [1,cos(polAng1),sin(polAng1),0;...
    cos(polAng1),cos(polAng1)^2,cos(polAng1)*sin(polAng1),0;...
    sin(polAng1),cos(polAng1)*sin(polAng1),sin(polAng1)^2,0;...
    0,0,0,0]/2;
P2 = [1,cos(polAng2),sin(polAng2),0;...
    cos(polAng2),cos(polAng2)^2,cos(polAng2)*sin(polAng2),0;...
    sin(polAng2),cos(polAng2)*sin(polAng2),sin(polAng2)^2,0;...
    0,0,0,0]/2;

for index = 1:length(t)
arg1 = 2*(2*pi*f1.*t(index) + p1);
arg2 = 2*(2*pi*f2.*t(index) + p2);
Cr1 = cos(arg1);
Sr1 = sin(arg1);
Cr2 = cos(arg2);
Sr2 = sin(arg2);
Ret1 = [1,0,0,0;...
    0,Cr1.^2+Sr1.^2.*Cd1,Cr1.*Sr1.*(1-Cd1),-Sr1.*Sd1;...
    0,Cr1.*Sr1.*(1-Cd1),Sr1.^2+Cr1^2.*Cd1,Cr1.*Sd1;...
    0,Sr1.*Sd1,-Cr1.*Sd1,Sd1].';
Ret2 = [1,0,0,0;...
    0,Cr2.^2+Sr2.^2.*Cd2,Cr2.*Sr2.*(1-Cd2),-Sr2.*Sd2;...
    0,Cr2.*Sr2.*(1-Cd2),Sr2.^2+Cr2^2.*Cd2,Cr2.*Sd2;...
    0,Sr2.*Sd2,-Cr2.*Sd2,Sd2].';
I(index) = [1,0,0,0]*P2*Ret2*M*Ret1*P1*[1;0;0;0];
end