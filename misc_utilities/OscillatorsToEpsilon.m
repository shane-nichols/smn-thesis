function epsilon = OscillatorsToEpsilon(Amp,E,Gamma,Theta,Phi,Lam)

cphi = cos(Phi);
sphi = sin(Phi);
cthet = cos(Theta);
sthet = sin(Theta);
eV = (1239.8./Lam);
epsilon = zeros(length(Lam),3,3);
for index = 1:length(eV)
    Osc = Amp./(E.^2 - eV(index).^2 + eV(index).*Gamma*1i);
    temp = zeros(3);
    for n = 1:length(Osc)
        Vect = [cthet(n).*sphi(n),sthet(n).*sphi(n),cphi(n)];
        temp = kron(Vect,Vect')*Osc(n)+temp;
    end
    epsilon(index,:,:) = temp;
end
end

