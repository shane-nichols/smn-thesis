function p = sl2cInvert(G)
% Caution: this assumes that det(G) = 1. If this is not the case, use
% gl2cInvert instead. 
K = 1/sqrt(det(G));
T = acos(K*(G(1,1) + G(2,2))/2);
O = T*K/sin(T);
%O = T/sqrt(1-((G(1,1) + G(2,2))/2).^2); <- may be better for symbolic stuff
p(1) = 1i*O*(G(1,1) - G(2,2));
p(2) = 1i*O*(G(1,2)+G(2,1));
p(3) = O*(G(2,1) - G(1,2));

end