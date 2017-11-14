function p = gl2cInvert(G)
% returns the principal element of the Lie algebra sl(2,C) corresponding to
% the element 'G' of the Lie group SL(2,C). By principal I mean first
% branch of the matrix log.
K = 1/sqrt(det(G));
T = acos(K*(G(1,1) + G(2,2))/2);
O = T*K/sin(T);
p(1) = 2*1i*log(1/K);
p(2) = 1i*O*(G(1,1) - G(2,2));
p(3) = 1i*O*(G(1,2)+G(2,1));
p(4) = O*(G(2,1) - G(1,2));

end