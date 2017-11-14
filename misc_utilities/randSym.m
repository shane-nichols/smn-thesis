function m = randSym(n)
% randon symmetric matrix of size n x n
m = rand(n);
m=triu(m,1)+triu(m,0).' ;
end

