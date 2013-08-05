clear all
close all
clc

mex /Users/sebastien/Desktop/WTGen/SplineGen/SplineCMex.c
!cp /Users/sebastien/Desktop/WTGen/SplineGen/SplineCMex.c /Users/sebastien/Desktop/WTGen/TestHessianApprox

R = 60;
W = 10;

[beta,lambda] = meshgrid(-1.9:0.1:20,2.1:0.1:15);
[m,n] = size(beta);

Mem = [];
for i = 1:m
    for j = 1:n
        beta_ij = beta(i,j);
        lambda_ij = lambda(i,j);

        X = SplineCMex(beta_ij,lambda_ij);

        Cp               = X(1);
        dCpdbeta         = X(2);
        dCpdlambda       = X(3);
        d2Cpdbeta2       = X(4);
        d2Cpdlambda2     = X(5);
        d2Cpdlambdadbeta = X(6);

        CpHess = [d2Cpdbeta2       d2Cpdlambdadbeta;
                  d2Cpdlambdadbeta d2Cpdlambda2];

 
        CpTab(i,j) = X(1);
        if (Cp > 0)
            Mem = [Mem;beta_ij lambda_ij eig(-CpHess+0.2*eye(2)).' Cp];
        end
    end
end
[mineig,i] = min(Mem(:,3:4));
      
display('Minimal eigenvalue')
mineig
display('at beta, lambda = ')
Mem(i,1)
Mem(i,2)

figure(1);clf
mesh(beta,lambda,max(0,CpTab));hold on
plot3(Mem(i,1),Mem(i,2),Mem(i,5),'marker','o','markersize',10,'markerfacecolor','k','markeredgecolor','k')






