A = rand(1,4)*3; %random PEM amplitudes
p = rand(1,4)*pi; %random PEM phases
f = 1000*[50.231,43.313,47.431,61.421]; %set PEM frequencies
M = rand(4); %random Mueller matrix
times = linspace(0.01,1,1000);
for idx = 1:length(times)
t=linspace(0,times(idx),30000).'; %time points
I = PEMmakeI(f,p,A,M,t);
[Mout1,Mout2]=PEMdirectDemod(f,p,A,I,t);
er1(idx)=sum(sum((M-Mout1).^2))./16;
er2(idx)=sum(sum((M-Mout2).^2))./16;
end
plot(times,log10(er1),times,log10(er2))
h=axes;
semilogy(h,times,er1,times,er2);
set(h,'YMinorTick','on');


