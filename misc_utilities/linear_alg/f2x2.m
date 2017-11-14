function fA = f2x2(A,fun)
% Evaluates analytic functions on 2x2 matrices.
% A := any 2x2 real or complex matrix with non-degenerate eigenvalues. 
% If fun(x) with x = 0 is not defined, then A must be non-singular.
% fun := handle to a function. i.e., @exp, @sqrt, @(x) (x+5)/x^2
val1 = A(1,1) + A(2,2);
val2 = sqrt( (A(1,1) - A(2,2) )^2 + 4*A(1,2)*A(2,1) );
a1 = (val1 + val2)/2; % eigenvalues
a2 = (val1 - val2)/2;
fa1 = fun(a1);
fa2 = fun(a2);
fA = ( (a1*fa2- a2*fa1).*eye(2) + (fa1 - fa2).*A )./(a1 - a2);
end