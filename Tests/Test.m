clear all
close all
clc

mex /Users/sebastien/Desktop/WTGen/SplineGen/CpMex.c

R = 60;
W = 10;


beta = 15
omega = 6*2*pi/60;
lambda = R*omega/W

X = CpMex(beta,lambda);

Cp               = X(1)
dCpdbeta         = X(2);
dCpdlambda       = X(3);
d2Cpdbeta2       = X(4);
d2Cpdlambda2     = X(5);
d2Cpdlambdadbeta = X(6);

CpHess = [d2Cpdbeta2       d2Cpdlambdadbeta;
          d2Cpdlambdadbeta d2Cpdlambda2];

 
reg = 0.01      
-CpHess + reg*eye(2)      


eig(-CpHess + reg*eig(2))
      








