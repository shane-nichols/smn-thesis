function [freq,phase,coef] = PEMfreqBasis
syms a1 a2 a3 a4 freq phase coef

x=-5:2:5;
y=-4:2:4;
lx = length(x);
ly = length(y);

for i=1:lx
    for j=1:lx
        c1(i,j) = a1*x(i) + a2*x(j);
    end
end

for i=1:ly
    c2(i) = a1*y(i);
end

for i=1:lx
    for j=1:ly
        c3(i,j) = a1*x(i) + a2*y(j);
    end
end


for i=1:lx
    for j=1:lx
        c4(i,j) = a3*x(i) + a4*x(j);
    end
end

for i=1:lx
    for j=1:lx
        for k=1:lx
            for l=1:lx
                c5(i,j,k,l) = a1*x(i) + a2*x(j) + a3*x(k) + a4*x(l);
            end
        end
    end
end

for j=1:ly
    for k=1:lx
        for l=1:lx
            c6(j,k,l) = a1*y(j) + a3*x(k) + a4*x(l);
        end
    end
end

for i=1:lx
    for j=1:ly
        for k=1:lx
            for l=1:lx
                c7(i,j,k,l) = a1*x(i) + a2*y(j) + a3*x(k) + a4*x(l);
            end
        end
    end
end

for i=1:ly
    c8(i) = a4*y(i);
end

for j=1:lx
    for k=1:lx
        for l=1:ly
            c9(j,k,l) = a1*x(j) + a2*x(k) + a4*y(l);
        end
    end
end

for i=1:ly
    for j=1:ly
        c10(i,j) = a1*y(i) + a4*y(j);
    end
end

for j=1:lx
    for k=1:ly
        for l=1:ly
            c11(j,k,l) = a1*x(j) + a2*y(k) + a4*y(l);
        end
    end
end

for i=1:ly
    for j=1:lx
        c12(i,j) = a3*y(i) + a4*x(j);
    end
end

for i=1:lx
    for j=1:lx
        for k=1:ly
            for l=1:lx
                c13(i,j,k,l) = a1*x(i) + a2*x(j) + a3*y(k) + a4*x(l);
            end
        end
    end
end

for j=1:ly
    for k=1:ly
        for l=1:lx
            c14(j,k,l) = a1*y(j) + a3*y(k) + a4*x(l);
        end
    end
end

for i=1:lx
    for j=1:ly
        for k=1:ly
            for l=1:lx
                c15(i,j,k,l) = a1*x(i) + a2*y(j) + a3*y(k) + a4*x(l);
            end
        end
    end
end

freq = {0,c1(:),c2(:),c3(:),c4(:),c5(:),c6(:),c7(:),c8(:),c9(:),c10(:),c11(:),...
    c12(:),c13(:),c14(:),c15(:)};

clear c c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15

for i=1:lx
    for j=1:lx
        c1(i,j) = a1*x(i) + a2*x(j) + pi/2*( sign(x(i)) + sign(x(j)) );
    end
end

for i=1:ly
    c2(i) = a1*y(i);
end

for i=1:lx
    for j=1:ly
        c3(i,j) = a1*x(i) + a2*y(j) + pi/2*sign(x(i)) ;
    end
end


for i=1:lx
    for j=1:lx
        c4(i,j) = a3*x(i) + a4*x(j) + pi/2*( sign(x(i)) + sign(x(j)) );
    end
end

for i=1:lx
    for j=1:lx
        for k=1:lx
            for l=1:lx
                c5(i,j,k,l) = a1*x(i) + a2*x(j) + a3*x(k) + a4*x(l) + ...
                    pi/2*( sign(x(i)) + sign(x(j)) +sign(x(k)) +sign(x(l)));
            end
        end
    end
end

for j=1:ly
    for k=1:lx
        for l=1:lx
            c6(j,k,l) = a1*y(j) + a3*x(k) + a4*x(l) + pi/2*( sign(x(k)) + sign(x(l)) );
        end
    end
end

for i=1:lx
    for j=1:ly
        for k=1:lx
            for l=1:lx
                c7(i,j,k,l) = a1*x(i) + a2*y(j) + a3*x(k) + a4*x(l) + ...
                    pi/2*( sign(x(i)) + sign(x(k)) +sign(x(l)) );
            end
        end
    end
end

for i=1:ly
    c8(i) = a4*y(i);
end

for j=1:lx
    for k=1:lx
        for l=1:ly
            c9(j,k,l) = a1*x(j) + a2*x(k) + a4*y(l) + pi/2*( sign(x(j)) + sign(x(k)) );
        end
    end
end

for i=1:ly
    for j=1:ly
        c10(i,j) = a1*y(i) + a4*y(j);
    end
end

for j=1:lx
    for k=1:ly
        for l=1:ly
            c11(j,k,l) = a1*x(j) + a2*y(k) + a4*y(l) + pi/2*sign(x(j));
        end
    end
end

for i=1:ly
    for j=1:lx
        c12(i,j) = a3*y(i) + a4*x(j) + pi/2*sign(x(j));
    end
end

for i=1:lx
    for j=1:lx
        for k=1:ly
            for l=1:lx
                c13(i,j,k,l) = a1*x(i) + a2*x(j) + a3*y(k) + a4*x(l) + ...
                    pi/2*( sign(x(i)) + sign(x(j)) +sign(x(l)) );
            end
        end
    end
end

for j=1:ly
    for k=1:ly
        for l=1:lx
            c14(j,k,l) = a1*y(j) + a3*y(k) + a4*x(l) + pi/2*sign(x(l));
        end
    end
end

for i=1:lx
    for j=1:ly
        for k=1:ly
            for l=1:lx
                c15(i,j,k,l) = a1*x(i) + a2*y(j) + a3*y(k) + a4*x(l) ...
                    + pi/2*( sign(x(i)) + sign(x(l)) );
            end
        end
    end
end
phase = {0,c1(:),c2(:),c3(:),c4(:),c5(:),c6(:),c7(:),c8(:),c9(:),c10(:),c11(:),...
    c12(:),c13(:),c14(:),c15(:)};

% for i=1:16
% freq{i} = subs(c{i},[a1,a2,a3,a4],f);
% phase{i} = subs(c{i},[a1,a2,a3,a4],p);
% end


clear c c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15
x = abs(x);
y = abs(y);

for i=1:lx
    for j=1:lx
        c1(i,j) = besselj(x(i),a1)*besselj(x(j),a2);
    end
end

for i=1:ly
    c2(i) = besselj(y(i),a1);
end

for i=1:lx
    for j=1:ly
        c3(i,j) = besselj(x(i),a1)*besselj(y(j),a2);
    end
end


for i=1:lx
    for j=1:lx
        c4(i,j) = besselj(x(i),a3)*besselj(x(j),a4);
    end
end

for i=1:lx
    for j=1:lx
        for k=1:lx
            for l=1:lx
                c5(i,j,k,l) = besselj(x(i),a1)*besselj(x(j),a2)*besselj(x(k),a3)*besselj(x(l),a4);
            end
        end
    end
end

for j=1:ly
    for k=1:lx
        for l=1:lx
            c6(j,k,l) = besselj(y(j),a1)*besselj(x(k),a3)*besselj(x(l),a4);
        end
    end
end

for i=1:lx
    for j=1:ly
        for k=1:lx
            for l=1:lx
                c7(i,j,k,l) = besselj(x(i),a1)*besselj(y(j),a2)*besselj(x(k),a3)*besselj(x(l),a4);
            end
        end
    end
end

for i=1:ly
    c8(i) = besselj(y(i),a4);
end

for j=1:lx
    for k=1:lx
        for l=1:ly
            c9(j,k,l) = besselj(x(j),a1)*besselj(x(k),a2)*besselj(y(l),a4);
        end
    end
end

for i=1:ly
    for j=1:ly
        c10(i,j) = besselj(y(i),a1)*besselj(y(j),a4);
    end
end

for j=1:lx
    for k=1:ly
        for l=1:ly
            c11(j,k,l) = besselj(x(j),a1)*besselj(y(k),a2)*besselj(y(l),a4);
        end
    end
end

for i=1:ly
    for j=1:lx
        c12(i,j) = besselj(y(i),a3)*besselj(x(j),a4);
    end
end

for i=1:lx
    for j=1:lx
        for k=1:ly
            for l=1:lx
                c13(i,j,k,l) = besselj(x(i),a1)*besselj(x(j),a2)*besselj(y(k),a3)*besselj(x(l),a4);
            end
        end
    end
end

for j=1:ly
    for k=1:ly
        for l=1:lx
            c14(j,k,l) = besselj(y(j),a1)*besselj(y(k),a3)*besselj(x(l),a4);
        end
    end
end

for i=1:lx
    for j=1:ly
        for k=1:ly
            for l=1:lx
                c15(i,j,k,l) = besselj(x(i),a1)*besselj(y(j),a2)*besselj(y(k),a3)*besselj(x(l),a4);
            end
        end
    end
end

coef = {0,c1(:),c2(:),c3(:),c4(:),c5(:),c6(:),c7(:),c8(:),c9(:),c10(:),c11(:),...
    c12(:),c13(:),c14(:),c15(:)};
end
