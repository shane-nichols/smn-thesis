function I = PEMmakeI(f,p,A,M,t,noise)
% simulation of the 4-PEM light intensity, in vectorized form.
% transmission axis of input polizer is assumed to be along xy
% transmission axis of output polizer is assumed to be along yx
% PEMs are assumed to be at 0,45,45,0 deg.
% I = column vector of light intensity values at detector
% f = row vector of PEM frequencies of in Hz
% p = row vector of PEM phases in radians
% A = row vector of PEM modulation amplitudes in radians
% M = 4x4 test Mueller matrix
% t = time values to compute I in seconds
% noise = amount of noise to add with respect to signal level of 1. 

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
A = reshape(M.',16,1);
I = B*A./4;
I = I + noise*2*(rand(size(I))-0.5); % add the nosie
I = round(2^16*I) / (2^16); % this line applies the digitizer finite bit-depth

end