function dG = dSO3(p,x,n)
% n-th derivative d^nG/dx^n of a group element in SO3 formed by parameter
% vector p.
theta = norm(p);
p = p(:)./theta;
arg = theta*x + n*pi/2;
dG = theta.^n* ( cos(arg)*(eye(3) - p*p.') + sin(arg)*skew(p) );
end