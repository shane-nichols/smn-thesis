% this script symmetrizes the magnetoelectric pseudotensor under all
% crystallographic point groups


m = sym('m',[3,3]);

% 2 (2 || y)
% define the symmetry operators
ry = diag([-1 1 -1]);
sym_2 = (m + ry*m*ry.')/2

% m (m -| y)
% define the symmetry operators
ry = diag([1 -1 1]);
sym_m = (m - ry*m*ry.')/2  % inversion operators change the sign

% 222 
rx = diag([1 -1 -1]);
ry = diag([-1 1 -1]);
rz = diag([-1 -1 1]); 

sym1 = (m + rx*m*rx.')/2; % symmetrized under x2
sym2 = (sym1 + ry*sym1*ry.')/2; % symmetrized under y2
sym_222 = (sym2 + rz*sym2*rz.')/2 % symmetrizd under z2

% mm2
mx = diag([-1 1 1]);
my = diag([1 -1 1]);
rz = diag([-1 -1 1]);

sym1 = (m - mx*m*mx.')/2;
sym2 = (sym1 - my*sym1*my.')/2;
out_mm2 = (sym2 + rz*sym2*rz.')/2

% mmm
mx = diag([-1 1 1]);
my = diag([1 -1 1]);
mz = diag([-1 -1 1]);

sym1 = (m - mx*m*mx.')/2;
sym2 = (sym1 - my*sym1*my.')/2;
sym_mmm = (sym2 - mz*sym2*mz.')/2

% 4
r4z = [0 1 0 ; -1 0 0 ; 0 0 1];
sym_4 = (m + r4z*m*r4z.' + (r4z)^2*m*(r4z.')^2 + (r4z)^3*m*(r4z.')^3)./4

% -4
r4z = [0 1 0 ; -1 0 0 ; 0 0 -1];
sym_4bar = (m - r4z*m*r4z.' + (r4z)^2*m*(r4z.')^2 - (r4z)^3*m*(r4z.')^3)./4

% 422
sym1 = sym_4; % use previous result
sym2 = (sym1 + rx*sym1*rx)/2;
sym_422 = (sym2 + ry*sym2*ry)/2;

% 4mm
sym1 = sym_4; % use previous result
sym2 = (sym1 - mx*sym1*mx)/2;
sym_4mm = (sym2 - my*sym2*my)/2

% -42m (better thought of as -422)
sym1 = sym_4bar; % use previous result
sym2 = (sym1 + rx*sym1*rx)/2;
sym_4bar2m = (sym2 + ry*sym2*ry)/2

% 3
r3z = [-1/2,sqrt(3)/2,0 ; -sqrt(3)/2,-1/2,0 ; 0,0,1];
sym_3 = simplify((m + r3z*m*r3z.' + (r3z)^2*m*(r3z.')^2 )./3)

% 32 (better thought of as 322)
sym1 = sym_3;
sym2 = (sym1 + rx*sym1*rx)/2;
sym_322 = (sym2 + ry*sym2*ry)/2

%3m (better thought of as 3mm)
sym1 = sym_3;
sym2 = (sym1 - mx*sym1*mx)/2;
sym_3mm = (sym2 - my*sym2*my)/2

% 6
m = rand(3)
r6z = [1/2,sqrt(3)/2,0 ; -sqrt(3)/2,1/2,0 ; 0,0,1];
sym_6 = (m + r6z*m*r6z.' + (r6z)^2*m*(r6z.')^2 ...
           + (r6z)^3*m*(r6z.')^3 + (r6z)^4*m*(r6z.')^4 ...
           + (r6z)^5*m*(r6z.')^5)./6  % did this numerically because Matlab is dumb
m = sym('m',[3,3]);
sym_6 =[m(1,1) m(1,2) 0 ; -m(1,2) m(1,1) 0 ; 0 0 m(3,3)]; 


%622 same as 422 and 322
%6mm same as 4mm and 3m
