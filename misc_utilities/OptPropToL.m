function L = OptPropToL(LB,LD,LBp,LDp,CB,CD,k_iso)
k = - k_iso - sqrt(LD^2 + LDp^2 + CD^2);
L = [k -LD -LDp CD ; -LD k CB LBp ; -LDp -CB k -LB ; CD -LBp LB k];
end
