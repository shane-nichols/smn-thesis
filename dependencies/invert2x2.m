function M = invert2x2(M)  
% analytically inverts an array of 2x2 matrices in-place
% M is an N-dimensional array of 2x2 matrices of size 2,2,N,M,O...
% This is still faster than multinv( ) on my machine.
sz = size(M);
M = reshape(M, 4, []);
M([1,4],:) = M([4,1], :);
M([2,3],:) = -M([2,3],:);
M = M ./ (M(1,:).*M(4,:) - M(3,:).*M(2,:));
M = reshape(M, sz);
end