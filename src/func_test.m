%% 
% 码猿: 马叶椿
% 来源：原创
% 版本：v1.0-2022.09.14

%% 
function [test,Ns_test] = func_test(L,Dictionary_base,num_max,X,paramL)

coef = cell(1,num_max);
for j=1:num_max 
    coef{j} = OMP(Dictionary_base{j},X,paramL);
end

resignal = cell(1,num_max);
Pu_resingal = cell(1,num_max);
eigen_resignla = cell(1,num_max);
norm_eigen = zeros(1,num_max);
for j=1:num_max
    resignal{j} = Dictionary_base{j}*coef{j};
    R_resignal = resignal{j} * resignal{j}'/L;
    [U,eigen_resignla{j},~] = svd(R_resignal);
%     Pu_resingal{j} = U(:,1:j)*(U(:,1:j)'*U(:,1:j))^(-1)*U(:,1:j)';
    eigen_resignla{j} = diag(eigen_resignla{j});
    norm_eigen(j) = norm(eigen_resignla{j},1);
    Pu_resingal{j} = U*(U'*U)^(-1)*U';
end
[Ux,~,~]=svd(X*X'/L);
Pux = Ux*(Ux'*Ux)^(-1)*Ux';
difference = cell(1,num_max);
eigenvalue_difference = cell(1,num_max);
eigenvector_difference = cell(1,num_max);
error_eigen = zeros(1,num_max);
SRM = zeros(1,num_max);
for j=1:num_max
    difference{j} = resignal{j}-X;
    Rd = difference{j}*difference{j}'/L;
    [eigenvector_difference{j},eigenvalue_difference{j},~]=svd(Rd);
    eigenvalue_difference{j} = diag(eigenvalue_difference{j});
    error_eigen(j) = norm(eigenvalue_difference{j},1);
    Pud = eigenvector_difference{j}*(eigenvector_difference{j}'*eigenvector_difference{j})^(-1)*eigenvector_difference{j}';
    SRM(j) = trace(Pux*Pud);
end
% [test,Ns_test] = min(error_eigen);
[test,Ns_test] = max(norm_eigen);
end