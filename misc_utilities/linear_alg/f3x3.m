function fA = f3x3(A,fun)
% Evaluates analytic functions on 3x3 matrices.
% A := any 3x3 real or complex matrix with non-degenerate eigenvalues.
% If fun(x) with x = 0 is not defined, then A must be non-singular.
% fun := handle to a function. i.e., @exp, @sqrt, @(x) (x+5)/x^2
a = eig(A); % kindly ask Matlab for the eigenvalues (faster and easier).
fa1 = fun(a(1));
fa2 = fun(a(2));
fa3 = fun(a(3));
dif12 = (a(1)-a(2));
dif13 = (a(1)-a(3));
dif23 = (a(2)-a(3));
fdif12 = fa1 - fa2;
fdif23 = fa2 - fa3;
fdif13 = fa1 - fa3;

fA = ( -a(2)*a(3)*dif23*fa1 + a(1)*a(3)*dif13*fa2 - a(1)*a(2)*dif12*fa3 ).*eye(3) + ...
    (-a(3).^2*fdif12 + a(2).^2*fdif13 - a(1).^2*fdif23 ).* A + ...
    (a(3)*fdif12 + a(1)*fdif23 - a(2)*fdif13).* A^2;
fA = -fA./(dif12*dif13*dif23);
end