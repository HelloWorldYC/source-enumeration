%% 
function [MIC,Ns_MIC]=func_MIC(X,M,L_jk)
% X���۲�����
% M���۲����ݵ�ά�ȣ�����Ԫ����

MIC = zeros(1,M);

% �����ϵ������
EX = mean(X);
% VYY = (X-EX)*(X-EX)';
VYY = X*X'/L_jk;
V = (diag(VYY)).^(-1/2);
V = diag(V);
CYY = V*VYY*V;
[u,v] = eig(CYY);
v = diag(v);
% v = v + sqrt(sum(v));
v = sort(v,'ascend');
% R_1=cov(Y');
% V=(diag(R_1)).^(-0.5);
% V=diag(V);
% R=V*R_1*V;
    
for i=2:M
    MIC(i)=v(i)-v(i-1);
end

[MIC_max,Imic] = max(MIC);
Ns_MIC = M-Imic+1;

end