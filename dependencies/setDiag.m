function out = setDiag(diag)
% sets the diagonal in an array of matrices.
% diag is a 2D array where each row is the diagonal elements of a matrix
% out is a 3D array of diagonal matrices
out = zeros(size(diag, 1), size(diag, 1), size(diag, 2));
for i=1:size(diag, 1)
    out(i,i,:) = diag(i,:);
end
end