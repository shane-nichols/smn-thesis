function R = gyration2R(G)
% Converts a gyration tensor (i.e., magnetoelectric) to the optical
% rotation tensor on which the wavevector k operates.

% the equations below are equilvalent to this eigensystem approach, if that
% is more intuitive to you:
% [U,V] = eig(G);
% 
% R = zeros(3);
% R(1,1) = V(2,2) + V(3,3);
% R(2,2) = V(1,1) + V(3,3);
% R(3,3) = V(1,1) + V(2,2);
% 
% R = U*R/U;

R = -G;
R(1,1,:) = G(2,2,:) + G(3,3,:);
R(2,2,:) = G(1,1,:) + G(3,3,:);
R(3,3,:) = G(2,2,:) + G(1,1,:);
end