function psi = psiIso(layer, wavelengths, kx)

epsilon = materialLibIso(layer{1},wavelengths).^2;
psi = zeros(4,4,length(wavelengths));

%%% the vectorized code below is equivalent to this
% for j = 1:length(wavelengths)
%     q = sqrt(epsilon(j) - kx(j).^2);
%     t1 = q/epsilon(j);
%     t2 = 1/q;
%     psi(:,:,j) = [   t1,           0,          0,         -t2;...
%                       1,           0,          0,           1;...
%                       0,          t2,         -t1,          0;...
%                       0,           1,          1,           0];
% end

q = sqrt(epsilon - kx.^2);
t1 = q ./ epsilon;
t2 = 1 ./ q;
t3 = ones(1,length(wavelengths));

psi(1,1,:) = t1;
psi(1,4,:) = -t2;
psi(2,1,:) = t3;
psi(2,4,:) = t3;
psi(3,2,:) = t2;
psi(3,3,:) = -t1;
psi(4,2,:) = t3;
psi(4,3,:) = t3;
end
    
    