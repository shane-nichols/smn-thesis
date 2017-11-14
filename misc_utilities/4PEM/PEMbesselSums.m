function coef = PEMbesselSums(A)

coef = zeros(1,6);
s=0;
for j=1:2:5
    for k=1:2:5
        s = s + (besselj(j,A(1))*besselj(k,A(2))).^2;
    end
end
coef(1) = 4*s; %X1X2

s=0;
for j=2:2:6
    s = s + besselj(j,A(1)).^2;
end
coef(2) = besselj(0,A(1)).^2 + 2*s; %Y1

s=0;
for j=-5:2:5
    for k=-6:2:6
        s = s + (besselj(j,A(1))*besselj(k,A(2))).^2;
    end
end
coef(3) = s; %X1Y2


s=0;
for j=1:2:5
    for k=1:2:5
        s = s + (besselj(j,A(4))*besselj(k,A(3))).^2;
    end
end
coef(4) = 4*s; %X3X4

s=0;
for j=2:2:6
    s = s + besselj(j,A(4)).^2;
end
coef(5) = besselj(0,A(1)).^2 + 2*s; %Y4

s=0;
for j=-5:2:5
    for k=-6:2:6
        s = s + (besselj(j,A(4))*besselj(k,A(3))).^2;
    end
end
coef(6) = s; %X4Y3