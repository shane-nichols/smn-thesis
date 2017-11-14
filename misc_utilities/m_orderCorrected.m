function m = m_orderCorrected(M,p)
% this function returns the order corrected differential Mueller matrix 
% from the Mueller matrix M. 'p' is the order of anisotropy. This method
% uses the Jordan Conocial form as described in:

% Vincent Devlaminck and Razvigor Ossikovski, 
% Uniqueness of the differential Mueller matrix of uniform homogeneous media,
% Opt. Lett., 39, 2014, 3149--3152.

% this is less efficient than the analytical method that I developed. It
% gives the same answer more or less.

S = [1 0 0 0 ; 0 1 0 0; 0 0 1 1i ; 0 0 1i 1];
Sinv = inv(S);
deltaJ = zeros(4);
deltaJ(3,4) = -2*pi;
deltaJ(4,3) = 2*pi;

for i=1:size(M,3)
    [V,D] = eig(M(:,:,i));
    D = diag(D);
    %logm(S*diag(D)*Sinv) 
    Q = V*Sinv; %#ok<MINV>
    %t1 = sqrt(real(D(3)) + 2*imag(D(3)^2)); % ??
    t1 = log(norm(D(3)));
    %t2 = -atan(imag(D(3))/real(D(3)));
    t2 = angle(D(3));
    lnJ = [log(D(1)), 0, 0, 0; ...
           0, log(D(2)), 0, 0 ; ...
           0, 0,  t1 t2;...
           0, 0, -t2, t1];
    m(:,:,i) = Q*lnJ/Q + p.*(Q*deltaJ/Q);
end
           
       
    