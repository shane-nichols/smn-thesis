f1 = 1;
f2 = 1.25;
p1 = 3.234;
p2 = 1.234;
d1 = 1.5;
d2 = 1.6;
M = eye(4);
t = linspace(0,2,1001);
t = t(1:1000);
polAng1 = 1.4;
polAng2 = 2.123;

I = CRRmakeI2(f1,f2,p1,p2,d1,d2,M,t,polAng1,polAng2);
p = CRRcalibrate(f1,f2,I,t); % p = [p1,p2,d1,d2,polAng2];

QWP = [1,0,0,0;0,1,0,0;0,0,0,-1;0,0,1,0]; % fast axis along y ( LR_ang = 0 )
I2 = CRRmakeI2(f1,f2,p1,p2,d1,d2,QWP,t,polAng1,polAng2);
M2 = CRRharmonicDemod8(f1,f2,p(1),p(2),p(3),p(4),I2,t,0,p(5));
LR_ang = MMgetp(M2,'lbang');
p([1,2,5]) = p([1,2,5]) - LR_ang;

Mtest = rand(4);
I2 = CRRmakeI2(f1,f2,p1,p2,d1,d2,Mtest,t,polAng1,polAng2);
M2 = CRRharmonicDemod8(f1,f2,p1,p2,d1,d2,I2,t,polAng1,polAng2)