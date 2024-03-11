%% 利用稀疏重构前后的奇异值误差来做估计
% 码猿: 马叶椿
% 来源：原创
% 版本：v1.0-2022.06.26

%% 
function [MSRSE,Ns_MSRSE] = func_MSRSE(L,Dictionary_base,num_max,X,paramL)

coef = cell(1,num_max);
for j=1:num_max 
    coef{j} = OMP(Dictionary_base{j},X,paramL);
end

resignal = cell(1,num_max);
for j=1:num_max
    resignal{j} = Dictionary_base{j}*coef{j};
end

difference = cell(1,num_max);
eigenvalue_difference = cell(1,num_max);
error_eigen = zeros(1,num_max);
for j=1:num_max
    difference{j} = resignal{j}-X;
    Rd = difference{j}*difference{j}'/L;
    [~,eigenvalue_difference{j},~]=svd(Rd);
    eigenvalue_difference{j} = diag(eigenvalue_difference{j});
%     eigenvalue_difference{j} = sum(eigenvalue_difference{j})+eigenvalue_difference{j};
%         error_eigen(j) = sum(eigenvalue_error{j}.^2);
    error_eigen(j) = norm(eigenvalue_difference{j},1);
end
[MSRSE,Ns_MSRSE] = min(error_eigen);

end