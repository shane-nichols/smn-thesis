function M = invert4x4(M)
% numerically inverts an array of 4x4 matrices. I don't know of an
% analytical routine for 4x4 matrices that is faster than Matlab's inv()
% depriciated. Now using multinv( )
sz = size(M);
M = reshape(M, 4, 4, []);
for idx = 1:size(M, 3)
    M(:,:,idx) = inv(M(:,:,idx));
end
M = reshape(M, sz);
end