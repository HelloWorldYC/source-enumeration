%%
function [alpha,Ns_MSTDC]=func_MSTDC(X,M,L_jk)

STD = zeros(1,M);
ui = zeros(1,M);
alpha = zeros(1,M);

EX = mean(X);
% VYY = (X-EX)*(X-EX)';
VYY = X*X'/L_jk;
V = (diag(VYY)).^(-1/2);
V = diag(V);
CYY = V*VYY*V;
[u,v] = svd(CYY);
v = diag(v);
% v = v + sqrt(sum(v));
v = sort(v,'ascend');

for i=2:M
    ui(i)=(v(i)+v(i-1))/2;
    STD(i)=sqrt((v(i)-ui(i))^2+(v(i-1)-ui(i))^2);
end

for j=3:M
    alpha(j)=STD(j)-STD(j-1);
end

[MSTDC_max,Imstdc] = max(alpha);
Ns_MSTDC = M-Imstdc+1;

end