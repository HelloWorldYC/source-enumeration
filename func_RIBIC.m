%% 改进的BIC
% 码猿: 马叶椿
% 论文：On Accurate Source Enumeration: A New Bayesian Information Criterion（2021）
% 版本：v1.0-2022.02.27
%       v1.1-2022.02.28

%% 
function [RIBIC,Ns_RIBIC] = func_RIBIC(qita,m,n,R)
c = m/n;
miu = (1+sqrt(c))^2;
sigma = ((1+sqrt(c))^(4/3))/(n*sqrt(c));
delta = ((-(3/4)*log10(-8*pi*log10(1-qita)))^(2/3))*sigma+miu;
[U,S,V] = svd(R);
S = diag(S);
S = S + sqrt(sum(S));
fai = zeros(1,m-1);
tao = zeros(1,m-1);
for k=1:m-1
    fai(k) = ((m*n)+(m-n)*(k-1))*log10(1+((delta-1)/(m-k+1)))+(m-n)*log10(delta);
end

RIBIC = zeros(1,m-1);
for k=1:m-1
    tao(k) = (1/(m-k))*(sum(S(k+1:m)));
end
for k=1:m-1
    temp = 1;
    temp1 = 0;
    for j=1:k
       temp = temp*S(j);
       temp1 = temp1+fai(j)+m*log10(tao(k)/S(j));
    end
    RIBIC(k) = n*log10((tao(k)^(m-k))*temp)+temp1;
end

[RIBIC_min,Ns_RIBIC] = min(RIBIC);

end