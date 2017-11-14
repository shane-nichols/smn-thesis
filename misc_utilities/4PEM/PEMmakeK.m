function Int = PEMmakeK

xi = sym('xi');

A = rand(1,4)*3; %random PEM amplitudes
freq = 1000*[50.231,43.313,47.431,61.421]; %set PEM frequencies
T=1;

clear B 
B = @(xi) X_f(xi,A(1),T).*X_f(xi,A(2),T);
Int = int(B,xi,0,400000);

% B(:,3) = @(xi) Y_f(xi,A(1),T);
% B(:,4) = @(xi)X_f(xi,A(1),T).*Y(xi,A(2),T);
% B(:,5) = @(xi) -X_f(xi,A(3),T).*X_f(xiA(4),T);
% B(:,9) = @(xi) -Y_f(xi,A(4),T);
% B(:,13) = @(xi) X_(xi,(4),T).*Y_f(xi,(3),T);
% B(:,6) = B(:,2).*B(:,5);
% B(:,7) = B(:,3).*B(:,5);
% B(:,8) = B(:,4).*B(:,5);
% B(:,10) = B(:,2).*B(:,9);
% B(:,11) = B(:,3).*B(:,9);
% B(:,12) = B(:,4).*B(:,9);
% B(:,14) = B(:,2).*B(:,13);
% B(:,15) = B(:,3).*B(:,13);
% B(:,16) = B(:,4).*B(:,13);
% M = K*(B.')*I./length(t);


    function X = X_f(xi,A,T)
        X = 0;
        for idx = 1:2:7;
            arg =  2*pi*T*(xi-freq(1)*idx) ;
            X = X + 2 * besselj(idx,A) * sin(arg)/arg;
        end
    end

    function Y = Y_f(xi,A,T)
        arg = 2*pi*T*xi;
        Y = J(0,A)*sin(arg)/arg;
        for idx = 2:2:8;
            arg =  2*pi*T*(xi-freq(1)*idx) ;
            Y = Y + 2 * besselj(idx,A) * sin(arg)/arg;
        end
    end


end