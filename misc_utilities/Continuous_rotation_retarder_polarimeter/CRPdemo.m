
% f1 = rotation frequency of the polarizer, in Hz
% f2 = rotation frequency of the analyzer, in Hz
% p1 = phase (initial angle) of the polarizer, in rad
% p2 = phase (initial angle) of the analyzer, in rad
% r1 = degree of linear polarization in the incident light (polarization
%      bias). Can range from 0 to 1.
% phi1 = angle of the polarization ellipse on the incident light, in rad.
% r2 = degree of linear polarization bias in the detection system.
%      Can range from 0 to 1.
% phi2 = angle of the linear polarization bias in detection system, in rad.
% M = 4x4 (or 3x3 sub-block) Mueller matrix of the sample.
% t = array of measurement times

% I = waveform of intensity values at detector

f1=1;
f2=1.25;
p1=rand;
p2=rand;
r1=0.1;
r2=0.13;
phi1=rand;
phi2=rand;
M=randn(4);
t=linspace(0,2,1001);
t=t(1:1000);


% make the waveform
I = CRPmakeI(f1,f2,p1,p2,r1,r2,phi1,phi2,M,t);
% demodulation the Mueller matrix elements from the waveform
CRPdemod(f1,f2,p1,p2,r1,r2,phi1,phi2,I,t)

M = ones(4);
noise = randn(size(t))/100;
ret = 80:160;
er = zeros(4,4,length(ret));
cAr = zeros(1,length(ret));
idx = 1;
for i=ret
d1 = i*pi/180;
[I,~,c] = CRRmakeI(f1,f2,p1,p2,d1,d1,ones(4),t);
cAr(idx) = c;
er(:,:,idx) = sqrt((M - CRRlsDemodNew(f1,f2,p1,p2,d1+0.01,d1+0.01,I,t)).^2);
idx = idx + 1;
end
plotter(ret,cAr);
MMplot(ret,er);

