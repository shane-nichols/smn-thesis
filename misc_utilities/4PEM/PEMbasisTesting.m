
syms a1 a2 a3 a4
[freqFun,phaseFun,coefFun] = PEMfreqBasis;

for i=1:16
freq{i} = double(subs(freqFun{i},[a1,a2,a3,a4],f));
phase{i} = double(subs(phaseFun{i},[a1,a2,a3,a4],p));
end

for i=1:16
coef{i} = double(subs(coefFun{i},[a1,a2,a3,a4],A));
end
coef{1} = 1;

B1 = zeros(16,length(t));
B2 = zeros(16,length(t));
for j=1:16
    for i=1:length(t)
        B2(j,i) = sum(coef{j}.*cos(2*pi*freq{j}*t(i) + phase{j}));
        B1(j,i) = sum(coef{j}.*mag(exp(1i*2*pi*freq{j}*t(i))) );
    end
end

% for j=1:16
%     for i=1:length(t)
%         % B2(j,i) = sum(coef{j}.*cos(2*pi*freq{j}*t(i) + phase{j} + pi/2));
%         B2(j,i) = sum(coef{j}.*cos(2*pi*freq{j}*t(i) + pi/2));
%     end
% end
% B2(1,:) = ones(1,length(t));
%B1 = sqrt(real(B1).^2 + imag(B1).^2);
Mout = (reshape(inv(B1*B1.')*B1*I,[4,4]).')*4;
Mout2 = sqrt(real(Mout).^2 + imag(Mout).^2);