function [p_out,phase_out] = PEMphaseVoltCali(t,f,p)
% p_out = [m,b,s] array of fitting values.
% phase_out = phase of the PEM
% this function demostrates how to find a linear relation relating the
% PEM voltage to the amplitude of modulation. 
volts = 0:0.01:2; % create an array of voltages to apply to the PEM
Amps = 0.045 + 2.1*volts; % convert volts to amps using a linear equation.
                            % the values b = 0.045 and m = 2.1 are what we are
                            % trying to find.
                            
for i = 1:length(Amps)
I = 100*(1 + 0.95*sin( Amps(i)*sin(2*pi*t*f(1)+p(1)) )); % simulaiton of the waveform with scale factor
% c = 100;
C1 = sum(exp( 1i*2*pi*t*f(1)).*I)./length(t); % get amplitude of C_v1
C2 = sum(exp( 1i*2*pi*t*3*f(1)).*I)./length(t); % get amplitude of C_v2
[phase1(i),mag1(i)] = cart2pol(real(C1),imag(C1)); % convert to mag and phase
[phase2(i),mag2(i)] = cart2pol(real(C2),imag(C2)); % convert to mag and phase
end

% add pi to any phases less than zero, average over the phases, then subtract from 
% pi/2.
phase_out = pi/2 - sum(phase1+pi*(phase1<0))./length(phase1)
                                         

figure
plot(volts,mag1,volts,mag2) % plot the magnitudes

p0 = [1,0,100];  % define initial parameters vector with slope of 1 and offset of 0
% and scale factor 100 that one can estimate by looking at the plotted data.

% perform non-linear least-square regression to determine parameters. 
p_out = lsqcurvefit(@(param,volts)fitMags(param,volts),p0,volts,[mag1;mag2]);
end

function mags_out = fitMags(param,volts) % model function
    Amps = volts*param(1)+param(2); % convert volts to amps
    mags_out = param(3)*abs([besselj(1,Amps) ; besselj(3,Amps)]); %array to compare to data
end
