function M = skew(n)
% makes a skew matrix from the vector n
% generator of SO(3)
M = [0,-n(3),n(2);n(3),0,-n(1);-n(2),n(1),0];
end