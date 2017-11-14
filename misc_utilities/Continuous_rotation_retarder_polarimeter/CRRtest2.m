function CRRtest2(M)
polAng0 = 0;
d0=1.5;
p0=0;
p1=0;
f0=1;
f1=1.25;
t = linspace(0,2-2/1000,1000);
alpha = 2*(2*pi*f0.*t + p0);
beta = 2*(2*pi*f1.*t + p1);
Cd0 = cos(d0);
Sd0 = sin(d0);
for n=1:1000
    polAng0 = n*2*pi/1000; 
    Xdata(n)=polAng0;
    P0 = [1,cos(polAng0),sin(polAng0),0;...
        cos(polAng0),cos(polAng0)^2,cos(polAng0)*sin(polAng0),0;...
        sin(polAng0),cos(polAng0)*sin(polAng0),sin(polAng0)^2,0;...
        0,0,0,0];
for index=1:length(t)
arg0 = 2*(2*pi*f0.*t(index) + p0);
arg1 = 2*(2*pi*f1.*t(index) + p1);
Cr0 = cos(arg0);
Sr0 = sin(arg0);
Cr1 = cos(arg1);
Sr1 = sin(arg1);
Ret0 = [1,0,0,0;...
    0,Cr0.^2+Sr0.^2.*Cd0,Cr0.*Sr0.*(1-Cd0),-Sr0.*Sd0;...
    0,Cr0.*Sr0.*(1-Cd0),Sr0.^2+Cr0^2.*Cd0,Cr0.*Sd0;...
    0,Sr0.*Sd0,-Cr0.*Sd0,Sd0];
Mout = M*Ret0*P0;
I(index) = Mout(1,1);
end
B(1,n) = sum(I)/2;
B(2,n) = sum(I.*exp(1i*alpha));
B(3,n) = sum(I.*exp(1i*2*alpha));
end
Br = real(B);
Bi = imag(B);

for hand=1:3
    plot_h(hand)=subplot(3,1,hand);
    hold on
end

for j = 1:3
        plot(plot_h(j),Xdata,squeeze(Br(j,:)),Xdata,squeeze(Bi(j,:)))
end
axis(plot_h,'tight')