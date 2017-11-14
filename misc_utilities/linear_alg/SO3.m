function G = SO3(p)
p0 = norm(p);
if p0==0
    G=eye(3);
else
p = p(:)./p0;
G = cos(p0)*eye(3) + sin(p0)*skew(p) + (1 - cos(p0))*(p*p.');
end
end