function [M,M2,K,K2] = PEMdirectDemod2(f,p,A,I,t)

arg = 2*pi*t*f; %arguments without phases
for idx = 1:length(t)
    arg(idx,:) =  arg(idx,:)+p; %add the phases
end
arg = sin(arg); %take sine
for idx = 1:length(t)
    arg(idx,:) =  arg(idx,:).*A; %multiply by amplitudes
end
X = sin(arg);
Y = cos(arg);

sXA1 = (1-besselj(0,2*A(1)))./2;
sYA1 = (1+besselj(0,2*A(1)))./2;
sXA2 = (1-besselj(0,2*A(2)))./2;
sYA2 = (1+besselj(0,2*A(2)))./2;
sXA3 = (1-besselj(0,2*A(3)))./2;
sYA3 = (1+besselj(0,2*A(3)))./2;
sXA4 = (1-besselj(0,2*A(4)))./2;
sYA4 = (1+besselj(0,2*A(4)))./2;
j04 = besselj(0,A(4));
j01 = besselj(0,A(1));
temp = sYA1 - j01^2;
K_G = 2*[sYA1./temp,0,-j01/temp,0;0,1/(sXA1*sXA2),0,0;-j01/temp,0,1/temp,0;0,0,0,1/(sXA1*sYA2)];
temp = sYA4 - j04^2;
K_A = 2*[sYA4./temp,0,j04/temp,0;0,1/(sXA4*sXA3),0,0;j04/temp,0,1/temp,0;0,0,0,1/(sXA4*sYA3)];

K =  kron(K_A,K_G)./length(t);

B_G(:,1) = ones(1,length(t)); 
B_G(:,2) = X(:,1).*X(:,2);
B_G(:,3) = Y(:,1);
B_G(:,4) = X(:,1).*Y(:,2);
B_A(:,1) = ones(1,length(t)); 
B_A(:,2) = -X(:,3).*X(:,4);
B_A(:,3) = -Y(:,4);
B_A(:,4) = X(:,4).*Y(:,3);
B_A = B_A./2;
B_G = B_G./2;
% K_G = inv(B_G.'*B_G)*sqrt(length(t));
% K_A = inv(B_A.'*B_A)*sqrt(length(t));

K2 = kron(K_A,K_G);

M = zeros(4);
for i=1:length(t)
M = M + I(i) .* (B_G(i,:).'*B_A(i,:));
end
M = (K_G*M*K_A).'.*4./length(t);

B = zeros(length(t),16);
B(:,1) = ones(1,length(t)); 
B(:,2) = X(:,1).*X(:,2);
B(:,3) = Y(:,1);
B(:,4) = X(:,1).*Y(:,2);
B(:,5) = -X(:,3).*X(:,4);
B(:,9) = -Y(:,4);
B(:,13) = X(:,4).*Y(:,3);
B(:,6) = B(:,2).*B(:,5);
B(:,7) = B(:,3).*B(:,5);
B(:,8) = B(:,4).*B(:,5);
B(:,10) = B(:,2).*B(:,9);
B(:,11) = B(:,3).*B(:,9);
B(:,12) = B(:,4).*B(:,9);
B(:,14) = B(:,2).*B(:,13);
B(:,15) = B(:,3).*B(:,13);
B(:,16) = B(:,4).*B(:,13);
M2 = K*(B.')*I;
M2 = reshape(M2,4,4);

end