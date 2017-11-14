function M = Mscatter(u, phi)
% makes a Mueller matrix for elastic scattering of a dipole along the
% direction of the 3-vector u. phi is the scattering angle.
c = u(1)*cos(phi) + u(3)*sin(phi);
J = kron([c;u(2)],[u(1),u(2)]);
A = [1,0,0,1;1,0,0,-1;0,1,1,0;0,1i,-1i,0]/sqrt(2);
M = A*kron(J,conj(J))*A';
end