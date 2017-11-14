function J = OptPropToJ(LB,LD,LBp,LDp,CB,CD,n,k)
X = n-1i*k;
L = LB-1i*LD;
Lp = LBp-1i*LDp;
C = CB-1i*CD;
% N = (-1i/2).*[X+L,Lp+1i*C;Lp-1i*C,X-L];
% compute J = expm(N) analytically;
T = sqrt(L.^2 + Lp.^2 + C.^2);
cT = cos(T/2);
if T == 0
    sT = 0;
else
    sT = sin(T/2)./T;
end
LsT = 1i*L.*sT;
LpsT = 1i*Lp.*sT;
CsT = C.*sT;
J(1,1,:) = cT - LsT;
J(2,2,:) = cT + LsT;
J(1,2,:) = CsT - LpsT;
J(2,1,:) = - CsT - LpsT;
J = exp(-1i*X/2).*J;
end