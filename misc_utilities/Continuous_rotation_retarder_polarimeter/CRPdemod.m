function M = CRPdemod(f1,f2,p1,p2,r1,r2,phi1,phi2,I,t)

B = zeros(length(t),9);

for index = 1:length(t)
    
    arg1 = 2*(2*pi*f1.*t(index) + p1);
    arg2 = 2*(2*pi*f2.*t(index) + p2);
    
    P1 = [(r1*cos(2*phi1 - 2*arg1))/2 + 1/2 ;...
        cos(2*arg1)/2 + (r1*cos(2*phi1 - 4*arg1))/4 + (r1*cos(2*phi1))/4 ;...
        sin(2*arg1)/2 + (r1*sin(2*phi1))/4 - (r1*sin(2*phi1 - 4*arg1))/4 ];
    
    P2 = [(r2*cos(2*phi2 - 2*arg2))/2 + 1/2 ,...
        cos(2*arg2)/2 + (r2*cos(2*phi2 - 4*arg2))/4 + (r2*cos(2*phi2))/4 ,...
        sin(2*arg2)/2 + (r2*sin(2*phi2))/4 - (r2*sin(2*phi2 - 4*arg2))/4 ];
    
    B(index,:) = kron(P1.',P2);
end
M = reshape(pinv(B)*I,3,3);