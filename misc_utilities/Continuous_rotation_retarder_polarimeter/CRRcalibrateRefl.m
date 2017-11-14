function [p0,p1,d0] = CRRcalibrateRefl(f0,f1,d2,I,t)
% Proof of Equation 37. on line 35.
I=I';
alpha = 2*(2*pi*f0.*t);
beta = 2*(2*pi*f1.*t);

C(1) = sum(I)/2;
C(2) = sum(I.*exp(1i*alpha));
C(3) = sum(I.*exp(1i*beta));
C(4) = sum(I.*exp(2*1i*alpha));
C(5) = sum(I.*exp(2*1i*beta));
C(6) = sum(I.*exp(1i*(alpha+beta)));
C(7) = sum(I.*exp(1i*(alpha-beta)));
C(8) = sum(I.*exp(1i*(alpha-2*beta)));
C(9) = sum(I.*exp(1i*(2*alpha-beta)));
C(10) = sum(I.*exp(1i*2*(alpha+beta)));
C(11) = sum(I.*exp(1i*2*(alpha-beta)));

C = C/length(t)*2;
[theta, rho] = cart2pol(real(C),imag(C));
theta = mod(theta,2*pi);
theta*180/pi;
t1 = (pi-theta(6))/2;
t2 = -(theta(7)/2);
p1 = (t1-t2)/2*180/pi;
p0 = (t1+t2)/2*180/pi;

% p1 = (pi-theta(6)+theta(7))/4*180/pi;
% p0 = (pi-theta(6)-theta(7))/4*180/pi;
B=rho(11)+rho(10);
%acos((rho(6)-rho(4))/(rho(6)+rho(4))) %transmission
%acos((rho(5)-rho(4))/(rho(5)+rho(4))) %transmission
%acos((rho(6)-rho(3))/(rho(6)+rho(3))); %reflection
%acos((rho(5)-rho(3))/(rho(5)+rho(3))); %reflection
%val = -(rho(4)*(cos(1.6)-1)+B*(1+cos(1.6)))/(1-cos(1.6));
%acos( (rho(5)-B+val)/(rho(5)+B+val) );
%acos( (rho(4)-B+val)/(rho(4)+B+val) );
d0 = acos((rho(4)+rho(5)-2*B-(rho(4)+B)*cos(d2))/(rho(5)+B));
%d0 = (rho(5)-B)/(rho(5)+B);
%d1 = (rho(6)-B)/(rho(5)+B);

end

