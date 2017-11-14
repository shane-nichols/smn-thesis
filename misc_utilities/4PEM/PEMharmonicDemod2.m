function M = PEMharmonicDemod2(f,p,A,I,t)
% Advanced harmonic demodulation of a 4PEM waveform. A window function is
% applied and J0(A1) = J0(A4) != 0 is possible;
I = BlackmanNuttall(length(I)).'.*I; % apply window function
basis =...  % basis harmonics. 9 means not a member of the set.
[ 1 -1  9  9 ;...
  2  9  9  9 ;...
  1  0  9  9 ;...
  9  9  1 -1 ;...
  -1 -1  -1  -1 ;...
 -2  9  1  1 ;...
 1  0 -1  -1 ;...
  9  9  9  2 ;...
  1  1  9 -2 ;...
  2  9  9 -2 ;...
  1  0  9 -2 ;...
  9  9  0 -1 ;...
 -1  1  0 -1 ;...
  2  9  0 -1 ;...
 -1  0  0  -1 ];
bessel = ones(15,4);
for i = 1:15
    for j=1:4
        if basis(i,j) ~= 9; bessel(i,j) = besselj(abs(basis(i,j)),A(j)); end
    end
end
bessel = prod(bessel,2); 
basis(basis==9)=0;
freq = basis*f.'; 
phase = (basis*p.') - sum(rem(basis,2),2)*pi/2;
C = exp(1i*2*pi*freq*t.')*I./length(t); 
[phase2,mag] = cart2pol(real(C),imag(C));
DC = sum(I)./length(t); 
M = mag./(bessel);
M = reshape([DC ; M],[4,4]).';
% next 9 lines correct for J0(A1) = J0(A4) != 0
j01 = besselj(0,A(1));
j04 = besselj(0,A(4));
M(1,1) = M(1,1) + j01 * M(1,3) + j04 * M(3,1) - j04*j01 * M(3,3);
M(1,2) = M(1,2) - j04 * M(3,2);
M(1,3) = M(1,3) - j04 * M(3,3);
M(1,4) = M(1,4) - j04 * M(3,4);
M(2,1) = M(2,1) - j01 * M(2,3);
M(3,1) = M(3,1) - j01 * M(3,3);
M(4,1) = M(4,1) - j01 * M(4,3);
% apply the phases to get correct signs
M = 4*M.*reshape([1;sign(cos(phase + phase2))],[4,4]).';
% next lines corrects for the window's DC factor. It can be removed if the
% result will just be normalized by M(1,1)
M = M./0.3635819;
% apply sign changes for this instrument configuration
M([2,3],:) = -M([2,3],:); 
end