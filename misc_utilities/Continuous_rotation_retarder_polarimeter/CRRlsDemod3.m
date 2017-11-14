function [M,D] = CRRlsDemod3(f0,f1,p0,p1,d0,d1,I,t,Mr,Mt)
% Proof of Equation 34c. See line 31. 
B = zeros(16,length(t));
arg0 = 2*(2*pi*f0.*t + p0);
arg1 = 2*(2*pi*f1.*t + p1);
Cp0 = cos(arg0);
C2p0 = cos(2*arg0);
Sp0 = sin(arg0);
S2p0 = sin(2*arg0);
Cp1 = cos(arg1);
C2p1 = cos(2*arg1);
Sp1 = sin(arg1);
S2p1 = sin(2*arg1);
Cd0 = cos(d0);
Sd0 = sin(d0);
Cd1 = cos(d1);
Sd1 = sin(d1);
a = sin(d0);
b = sin(d1);
c = (1+cos(d0))/2;
d = -(1+cos(d1))/2;
e = (1-cos(d0))/2;
f = (cos(d1)-1)/2;
Dg=[0,0,0,0];
Da=[0;0;0;0];
D=zeros(16,1000);
Mat = kron(inv(Mt),transpose(inv(Mr)));
for index=1:length(t)
Dg(1)=(e-2*c*C2p0(index))/(2*e);
Dg(2)=C2p0(index)/e;
Dg(3)=S2p0(index)/e;
Dg(4)=Sp0(index)/a;
Da(1)=(f-2*d*C2p1(index))/(2*f);
Da(2)=C2p1(index)/f;
Da(3)=S2p1(index)/f;
Da(4)=Sp1(index)/b;
Da=Da*sqrt(2);
Dg=Dg*sqrt(2);
%Da=inv(Mt)*Da;
%Dg=Dg*inv(Mr);
%D(:,index) = Mat*kron(Da,transpose(Dg));
D(:,index) = kron(Da,transpose(Dg));
end
size(D)
M = reshape(D*I,4,4)'/length(t)*2;
end