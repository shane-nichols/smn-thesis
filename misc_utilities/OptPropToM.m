function M = OptPropToM(LB,LD,LBp,LDp,CB,CD,k_iso)
k = k_iso + sqrt(LD^2 + LDp^2 + CD^2);
J = OptPropToJ(LB,LD,LBp,LDp,CB,CD,0,k);
M = kron(J,conj(J));
Ainv = [0.5,0.5,0,0;0,0,0.5,-0.5*1i;0,0,0.5,0.5*1i;0.5,-0.5,0,0];
A = [1,0,0,1;1,0,0,-1;0,1,1,0;0,1i,-1i,0];
M = A*M*Ainv;
end
