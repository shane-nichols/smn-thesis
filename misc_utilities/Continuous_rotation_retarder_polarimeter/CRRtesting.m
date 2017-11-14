I = CRRmakeI(1,1.25,1.232,1.431,1.5435,1.6, Mt*Mref*Mr,t);

Mout = zeros(16,1000);
Mout2 = zeros(16,1000);
for index = 1:1000
    d0 = index*pi/4/1000-pi/8;
    p0 = 2*pi*index/1000;
    M = Mref*inv(Mr)*CRRlsDemod2(1,1.25,1.232,1.431,1.5435+d0,1.6,I,t);
    Mout(:,index) = reshape(M',1,16)./M(1,1);
    M = CRRlsDemod2(1,1.25,1.232,1.431,1.5435,1.6,I,t);
    Mout2(:,index) = reshape(M',1,16)./M(1,1);
    Xdata(index) = d0*180/pi^2;
end
figure
%Moutr = real(Mout2);
%Mouti = imag(Mout2);
for hand=1:16
    plot_h(hand)=subplot(4,4,hand);
    hold on
end

for j = 1:16
        plot(plot_h(j),Xdata,squeeze(Mout(j,:)),Xdata,squeeze(Mout2(j,:)))
end
axis(plot_h,'tight')
%adjust Ylimits
for index=1:length(plot_h)
    %ylim(plot_h(index),[-2,2]);
end
% 
% for hand=1:16
%     plot_h(hand)=subplot(4,4,hand);
%     hold on
% end
% 
% for j = 1:4
%     for l = 1:4
%         plot(plot_h(l+4*(j-1)),Xdata,squeeze(Mout(j,l,:)),Xdata,squeeze(Mout2(j,l,:)))
%     end
% end