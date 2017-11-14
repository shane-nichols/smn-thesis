function I = PEMmakeI2(f,p,A,M,t)
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
arg = 2*pi*t*f; %arguments without phases
for idx = 1:length(t)
    arg(idx,:) =  arg(idx,:)+p; %add the phases
end
arg = sin(arg); %take sine
for idx = 1:length(t)
    arg(idx,:) =  arg(idx,:).*A; %multiply by amplitudes
end
R = exp(1i*arg);
B = zeros(length(t),16);
B(:,1) = ones(1,length(t)); 
B(:,2) = 1i*R(:,1).*1i*R(:,2);
B(:,3) = R(:,1);
B(:,4) = 1i*R(:,1).*R(:,2);
B(:,5) = -1i*R(:,3).*1i*R(:,4);
B(:,9) = -R(:,4);
B(:,13) = 1i*R(:,4).*R(:,3);
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
I = B*A;
end