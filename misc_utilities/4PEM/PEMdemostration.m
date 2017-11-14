A = [2.4048 0.88 0.88 2.4048]; % general PEM amplitudes
p = rand(1,4)*pi; % random PEM phases
f = [42065.898, 49780.423, 59720.594, 47042.204]; %PEM frequencies
fs = 448000; % sampling rate
N = 1048576; % number of samples
t = transpose(0:(1/fs):((N-1)/fs)); % time array
M = [0.9758   -0.0088   -0.0077   -0.0134 ;...
    -0.0118   -0.3044    0.5257    0.7637 ;...
    -0.0095    0.7264    0.6345   -0.1470 ;...
     0.0092   -0.5757    0.5226   -0.5894 ]; % test Mueller matrix
M  = eye(4);
I = PEMmakeI(f,p,A,M,t,0);

[Mout1,Mout2] = PEMdirectDemod(f,p,A,I,t);
Mout3 = PEMharmonicDemod2(f,p,A,I,t);

% full inversion - analytic. Assumes waveform is infinite, kind of sucks.
Error1 = M - Mout1

% full inversion - numeric. This method is exact but requires accurate phases.

Error2 = M - Mout2

% inversion using 16 pure harmonics with window function applied.
% Phases can be known only approximately. 

Error3 = M - Mout3