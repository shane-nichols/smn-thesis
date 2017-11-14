function [M,M2] = PEMdirectDemod(f,p,A,I,t)

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
K_G = [sYA1./temp,0,-j01/temp,0;0,1/(sXA1*sXA2),0,0;-j01/temp,0,1/temp,0;0,0,0,1/(sXA1*sYA2)];
temp = sYA4 - j04^2;
K_A = [sYA4./temp,0,j04/temp,0;0,1/(sXA4*sXA3),0,0;j04/temp,0,1/temp,0;0,0,0,1/(sXA4*sYA3)];

K =  kron(K_A,K_G)./length(t); % analytic basis inversion matrix. 

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
M = K*(B.')*I; %analtyic inversion
M2 = reshape(inv(B.'*B)*B.'*I,[4,4]).'; % numeric Moore-Penrose Pseudoinverse
M = reshape(M,4,4).';
end