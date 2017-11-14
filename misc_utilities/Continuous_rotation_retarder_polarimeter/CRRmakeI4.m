function I = CRRmakeI4(f0,f1,p0,p1,M,t,pR1,pR2,pP1,pP2)

pP1iso = sqrt(imag(pP1).^2);
[~,P1] = gl2c(-1i*pP1iso,pP1);

pP2iso = sqrt(imag(pP2).^2);
[~,P2] = gl2c(-1i*pP2iso,pP2);

pR1iso = sqrt(imag(pR1).^2);
[~,R1] = gl2c(-1i*pR1iso,pR1);

pR2iso = sqrt(imag(pR2).^2);
[~,R2] = gl2c(-1i*pR2iso,pR2);

for index = 1:length(t)
arg0 = (2*pi*f0.*t(index) + p0);
arg1 = (2*pi*f1.*t(index) + p1);

Ret1 = MMrotate(R1,arg0);
Ret2 = MMrotate(R2,arg1);

I(index) = [1,0,0,0]*P2*Ret2*M*Ret1*P1*[1;0;0;0];
end
I=I';