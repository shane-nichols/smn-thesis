function [J,M] = sl2c(p)

T = sqrt(sum(p.^2));
%T = sqrt(L.^2 + Lp.^2 + C.^2);
cT = cos(T/2);
if T == 0;
    sT = 0;
else
    sT = sin(T/2)./T;
end
LsT = 1i*p(1).*sT;
LpsT = 1i*p(2).*sT;
CsT = p(3).*sT;
J(1,1,:) = cT - LsT;
J(2,2,:) = cT + LsT;
J(1,2,:) = - CsT - LpsT;
J(2,1,:) =  CsT - LpsT;
A = [1,0,0,1;1,0,0,-1;0,1,1,0;0,1i,-1i,0]./sqrt(2);
M = A*kron(J,conj(J))*A';

end