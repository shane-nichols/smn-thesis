function g = g_sl2c(p)
% Returns an element of the Lie algebra sl(2,C) where p is a parameter
% vector
g = -1i*[p(1) , p(2) - 1i*p(3) ; p(2) + 1i*p(3) , -p(1) ]/2;
end