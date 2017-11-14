function I = CRRmakeI3(f0,f1,p0,p1,d0,d1,M,t,polAng0,polAng1,polMag0,polMag1)
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


Cd0 = cos(d0);
Sd0 = sin(d0);
Cd1 = cos(d1);
Sd1 = sin(d1);

LD = cos(2*polAng0)*polMag0;
LDp = sin(2*polAng0)*polMag0;
G = -eye(4)*polMag0+SO31([1i*LD,1i*LDp,0]);
P0 = expm(G);

LD = cos(2*polAng1)*polMag1;
LDp = sin(2*polAng1)*polMag1;
G = -eye(4)*polMag1+SO31([1i*LD,1i*LDp,0]);
P1 = expm(G);

% polAng0 = polAng0*2;
% polAng1 = polAng1*2;
% P0 = [1,cos(polAng0),sin(polAng0),0;...
%     cos(polAng0),cos(polAng0)^2,cos(polAng0)*sin(polAng0),0;...
%     sin(polAng0),cos(polAng0)*sin(polAng0),sin(polAng0)^2,0;...
%     0,0,0,0]/2;
% P1 = [1,cos(polAng1),sin(polAng1),0;...
%     cos(polAng1),cos(polAng1)^2,cos(polAng1)*sin(polAng1),0;...
%     sin(polAng1),cos(polAng1)*sin(polAng1),sin(polAng1)^2,0;...
%     0,0,0,0]/2;

for index = 1:length(t)
arg0 = 2*(2*pi*f0.*t(index) + p0);
arg1 = 2*(2*pi*f1.*t(index) + p1);
Cr0 = cos(arg0);
Sr0 = sin(arg0);
Cr1 = cos(arg1);
Sr1 = sin(arg1);
Ret0 = [1,0,0,0;...
    0,Cr0.^2+Sr0.^2.*Cd0,Cr0.*Sr0.*(1-Cd0),-Sr0.*Sd0;...
    0,Cr0.*Sr0.*(1-Cd0),Sr0.^2+Cr0^2.*Cd0,Cr0.*Sd0;...
    0,Sr0.*Sd0,-Cr0.*Sd0,Sd0]
Ret1 = [1,0,0,0;...
    0,Cr1.^2+Sr1.^2.*Cd1,Cr1.*Sr1.*(1-Cd1),-Sr1.*Sd1;...
    0,Cr1.*Sr1.*(1-Cd1),Sr1.^2+Cr1^2.*Cd1,Cr1.*Sd1;...
    0,Sr1.*Sd1,-Cr1.*Sd1,Sd1];
I(index) = [1,0,0,0]*P1*Ret1*M*Ret0*P0*[1;0;0;0];
end
I=I';