function [I,P1out,P2out] = CRPmakeI(f1,f2,p1,p2,r1,r2,phi1,phi2,M,t)
% f1 = rotation frequency of the polarizer, in Hz
% f2 = rotation frequency of the analyzer, in Hz
% p1 = phase (initial angle) of the polarizer, in rad
% p2 = phase (initial angle) of the analyzer, in rad
% r1 = degree of linear polarization in the incident light (polarization
%      bias). Can range from 0 to 1.
% phi1 = angle of the polarization ellipse on the incident light, in rad.
% r2 = degree of linear polarization bias in the detection system.
%      Can range from 0 to 1.
% phi2 = angle of the linear polarization bias in detection system, in rad.
% M = 4x4 (or 3x3 sub-block) Mueller matrix of the sample.
% t = array of measurement times

% I = waveform of intensity values at detector

% P1 is the result of taking the first three elements of the vector formed
% by the matrix product:

    %  R(-theta)*M*R(theta)*S_in
    % S_in  = [ 1 ; r1*cos(2*phi1) ; r1*sin(2*phi2) ; 0]
    
% which is just the stokes vector leaving a rotatable polarization with
% incudent light having some degree of linear polarization. Note that the
% circular polarization on the incident light does not matter.

M = M(1:3,1:3);
I = zeros(length(t),1);
for index = 1:length(t)
    
    arg1 = 2*(2*pi*f1.*t(index) + p1);
    arg2 = 2*(2*pi*f2.*t(index) + p2);
    
    P1 = [(r1*cos(2*phi1 - 2*arg1))/2 + 1/2 ;...
        cos(2*arg1)/2 + (r1*cos(2*phi1 - 4*arg1))/4 + (r1*cos(2*phi1))/4 ;...
        sin(2*arg1)/2 + (r1*sin(2*phi1))/4 - (r1*sin(2*phi1 - 4*arg1))/4 ];
    
    P1out(:,index) = P1;
    
    P2 = [(r2*cos(2*phi2 - 2*arg2))/2 + 1/2 ,...
        cos(2*arg2)/2 + (r2*cos(2*phi2 - 4*arg2))/4 + (r2*cos(2*phi2))/4 ,...
        sin(2*arg2)/2 + (r2*sin(2*phi2))/4 - (r2*sin(2*phi2 - 4*arg2))/4 ];
    
    P2out(:,index) = P2;
    
    I(index) = P2*M*P1;
end
