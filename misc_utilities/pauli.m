function pauli = pauli(n)
% returns the nth Pauli matrix of the extended basis (includes identity) 

switch n
    case 1
        pauli = eye(2);
    case 2
        pauli = [1,0;0,-1];
    case 3
        pauli = [0,1;1,0];
    case 4
        pauli = [0,-1i;1i,0];
    otherwise
        error('n should be between 1 and 4')
end
end