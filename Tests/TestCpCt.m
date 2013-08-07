clear all
close all
clc

mex /Users/sebastien/Desktop/WTGen/SplineGen/CpMex.c
mex /Users/sebastien/Desktop/WTGen/SplineGen/CtMex.c


R = 60;
W = 10;
Cpmin_disp = 0;


[beta,lambda] = meshgrid(-1.9:0.1:25,2.1:0.1:15);
[m,n] = size(beta);

Mem = [];MemNeg = [];
for i = 1:m
    for j = 1:n
        beta_ij = beta(i,j);
        lambda_ij = lambda(i,j);

        X = CpMex(beta_ij,lambda_ij);
        Y = CtMex(beta_ij,lambda_ij);
        
   
        Cp               = X(1);
        dCpdbeta         = X(2);
        dCpdlambda       = X(3);
        d2Cpdbeta2       = X(4);
        d2Cpdlambda2     = X(5);
        d2Cpdlambdadbeta = X(6);

        CpTab(i,j) = X(1);
        CtTab(i,j) = Y(1);
        
   
    end
end
      

figure(1);clf
mesh(beta,lambda,max(Cpmin_disp,CpTab));hold on

figure(2);clf
mesh(beta,lambda,max(Cpmin_disp,CtTab));hold on





