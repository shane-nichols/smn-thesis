function coef = PEMbesselSums2(A)

coef = zeros(1,6);
% s=0;
% for j=1:2:20
%     for k=1:2:20
%         s = s + (besselj(j,A(1))*besselj(k,A(2))).^2;
%     end
% end
% coef(1) = 4*s; %X1X2

s=0;
for j=-21:2:21
    s = s + besselj(j,A(2)).^2;
end
coef(1) = s; %X1

s=0;
for j=-20:2:20
    s = s + besselj(j,A(1)).^2;
end
coef(2) = s; %Y1
coef(3) = coef(1)+coef(2);
s=0;
for j=-20:2:20
    for k=-20:2:20
        s = s + (besselj(j,A(1))*besselj(k,A(4)));
    end
end
coef(4)=s;
% s=0;
% for j=-21:2:21
%     for k=-20:2:20
%         s = s + (besselj(j,A(1))*besselj(k,A(2))).^2;
%     end
% end
% coef(3) = s; %X1Y2
end