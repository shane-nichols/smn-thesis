[~,Mtest] = gl2c(-1i*0.2,[1.2-0.3*1i,-0.2,0.1-0.1*1i]);
testWavelengths = 540:0.5:560; % width of 20 centerered at 550
retardances = interp1(Wavelength,Retardance0,testWavelengths);
t=linspace(0,4,2001);  % make time array from 0->4 seconds
t=t(1:2000); % remove the ducplicate point at 4 because it coincides with 0
I = CRRmakeI(1,1.25,0,0,Retardance0(4),Retardance0(4),Mtest,t); % simulate a waveform using retardances at 550
Mout = zeros(4,4,length(testWavelengths));
for i=1:length(testWavelengths)  % compute M from waveform using retardances at test wavelengths
    Mout(:,:,i) = CRRlsDemodNew(1,1.25,0,0,Retardance0(4),retardances(i),I,t);
end
Mavg = sum(Mout,3)./length(testWavelengths);
Mavg./Mavg(1,1) - Mtest./Mtest(1,1)


