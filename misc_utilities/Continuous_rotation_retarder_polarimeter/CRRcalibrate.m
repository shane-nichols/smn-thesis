function p = CRRcalibrate(f0,f1,I,t)
% Equation 26 - 27.
arg0 = 2*(2*pi*f0.*t);
arg1 = 2*(2*pi*f1.*t);
A(1) = sum(I.*exp(1i*(arg0+arg1)));
A(2) = sum(I.*exp(1i*(arg0-arg1)));
A(3) = sum(I.*exp(1i*2*(arg0+arg1)));
A(4) = sum(I.*exp(1i*2*(arg0-arg1)));
A(5) = sum(I.*exp(1i*2*(arg0)));
A(6) = sum(I.*exp(1i*2*(arg1)));
A = A/length(t)*2;
[theta, rho] = cart2pol(real(A),imag(A));
%theta*180/pi;
theta = mod(theta,2*pi);
polAng2 = ( theta(5)-theta(6)-theta(4) )/2;
%theta*180/pi;
p0 = (pi - theta(1) - theta(2))/4;
p1 = polAng2 + (theta(2) - theta(1) - pi)/4;
B=rho(3)-rho(4);
%acos((rho(6)-rho(4))/(rho(6)+rho(4))) %transmission
%acos((rho(5)-rho(4))/(rho(5)+rho(4))) %transmission
d0 = acos( 1 - ( (2*rho(4))/(rho(6)+rho(4)) )  );
d1 = acos(        2*rho(5) /(rho(5)+rho(4)) - 1);
%d0 = (rho(5)-B)/(rho(5)+B);
%d1 = (rho(6)-B)/(rho(5)+B);
p = [p0,p1,d0,d1,polAng2];

end

