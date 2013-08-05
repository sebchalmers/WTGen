clear all
close all
clc

syms Cp dCpdlambda dCpdbeta d2Cpdlambda2 d2Cpdbeta2 d2Cpdlambdadbeta real

syms W omega beta R real

Cost = -Cp*W^3;

lambda = R*omega/W;


Var = [beta;omega];

dlambda = jacobian(lambda, Var);
dbeta   = jacobian(beta,   Var);

dCp  = dCpdlambda*dlambda + dCpdbeta*dbeta;

dCost = jacobian(Cost,Var) + jacobian(Cost,Cp)*dCp;



d2Cost = jacobian(dCost,Var) + jacobian(dCost,Cp)*dCp + ...
         jacobian(dCost,dCpdlambda)*d2Cpdlambda2*dlambda + ...
         jacobian(dCost,dCpdlambda)*d2Cpdlambdadbeta*dbeta + ...
         jacobian(dCost,dCpdbeta)*d2Cpdlambdadbeta*dlambda + ...
         jacobian(dCost,dCpdbeta)*d2Cpdbeta2*dbeta;
     

List_of_variables = {'dCost', 'd2Cost'};
File_name = 'Hessian.m'

index_equ = 1;
for k = 1:length(List_of_variables) 
Var_name = List_of_variables{k};
size_var = size(eval(Var_name));
for i = 1:size_var(1)
    for j = 1:size_var(2)
        Var_name_ij = strcat(Var_name,'(',num2str(i),',',num2str(j),')');
        Equation = strcat(Var_name_ij,' = ',char(eval(Var_name_ij)),';');
        string2eval = strcat('Equ',num2str(index_equ),' = Equation;');
        eval(string2eval);
        index_equ = index_equ + 1;
    end
end
end
fid = fopen([File_name],'wt');
for k = 1:index_equ-1
    eval(strcat('fprintf(fid,''%s \n'',Equ',num2str(k),');'));
end
fclose(fid);
clear Equ*   
