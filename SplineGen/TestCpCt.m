clear all
close all
clc

mex /Users/sebastien/Desktop/WTGen/SplineGen/CpMex.c
mex /Users/sebastien/Desktop/WTGen/SplineGen/CtMex.c

run MLSModel

R = 60;
W = 10;
Cpmin_disp = 0;


[beta,lambda] = meshgrid(-4.9:0.1:25,2.1:0.1:15);
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

        CpSplined(i,j) = X(1);
        CtSplined(i,j) = Y(1);
        
   
    end
end
      

figure(1);clf
mesh(beta,lambda,max(Cpmin_disp,CpSplined));hold on
mesh(Cptable(1,2:end),Cptable(2:end,1),max(Cpmin_disp,Cptable(2:end,2:end)));hold on

figure(2);clf
mesh(beta,lambda,CtSplined);hold on
mesh(Cptable(1,2:end),Cptable(2:end,1),Cttable(2:end,2:end));hold on





