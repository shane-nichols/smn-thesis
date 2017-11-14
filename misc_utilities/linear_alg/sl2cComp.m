function v = sl2cComp(b,a)
na = sqrt(sum(a.^2));
ea = a(:)./na;
nb = sqrt(sum(b.^2));
eb = b(:)./nb;
na=na/2;
nb=nb/2;
sa = sin(na);
ca = cos(na);
cb = cos(nb);
sb = sin(nb);
nv = 2*acos(ca*cb - sa*sb*(ea.'*eb));
v = sa*cb*ea + ca*sb*eb + sa*sb*cross(eb,ea);
v = (nv./sqrt(sum(v.^2))).*v;

% simplify(sa*cb*ea + ca*sb*eb + sa*sb*cross(eb,ea))
% simplify(ca*cb - sa*sb*(ea.'*eb))
end