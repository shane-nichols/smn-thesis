function out = bigKron(a)
% Perfors kron(m, conj(m)) for array of matrices
% a is an N-dimensional array of 2x2 matrices of size 2,2,N,M,O...
% out is an array of size 4,4,N,N,O...
sz = size(a);
a = reshape(a, 2, 2, []);
out = zeros(4, 4, size(a, 3));
out(1,1,:) = a(1,1,:) .* conj(a(1,1,:));
out(1,2,:) = a(1,1,:) .* conj(a(1,2,:));
out(1,3,:) = conj( out(1,2,:) );
out(1,4,:) = a(1,2,:) .* conj(a(1,2,:));
out(2,1,:) = a(1,1,:) .* conj(a(2,1,:));
out(2,2,:) = a(1,1,:) .* conj(a(2,2,:));
out(2,3,:) = a(1,2,:) .* conj(a(2,1,:));
out(2,4,:) = a(1,2,:) .* conj(a(2,2,:));
out(3,1,:) = conj( out(2,1,:) );
out(3,2,:) = conj( out(2,3,:) );
out(3,3,:) = conj( out(2,2,:) );
out(3,4,:) = conj( out(2,4,:) );
out(4,1,:) = a(2,1,:) .* conj(a(2,1,:));
out(4,2,:) = a(2,1,:) .* conj(a(2,2,:));
out(4,3,:) = conj( out(4,2,:) );
out(4,4,:) = a(2,2,:) .* conj(a(2,2,:));
out = reshape(out, [4, 4, sz(3:end)] );
end