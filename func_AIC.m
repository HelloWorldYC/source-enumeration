%%
function  [AIC,Ns_AIC] = func_AIC(M,L,D)
% AIC函数的实现
AIC=zeros(1,M-1);
for n=1:M-1
    temp0=sum(D(n+1:M))/(M-n);
    temp1=1;
    for ii=n+1:M
        temp1=temp1*D(ii);
    end
    temp1=temp1^(1/(M-n));
    AIta=temp0/temp1;   
    AIC(n)=2*L*(M-n)*log10(AIta)+2*n*(2*M-n);
end
[AIC_min,Ns_AIC] = min(AIC);
% Ns_AIC = M-(Ns_AIC-1)+1;
% Ns_AIC = Ns_AIC-1;
end

%%
% function [num]=func_AIC(R,N,M,D)
% %AIC准则估计信源数目
% %input R 阵元接收信号的协方差矩阵
% %input N 采样数
% %input M 阵元数
% %input D 1：一维线阵 2：二维方阵
% %output num 信源数目
% if D==1
%     [u,s,v] = svd(R);
%     sd = diag(s);
%     sd=sort(sd,'descend');
%     a = zeros(1,M);
%     for m = 0:M-1
%         negv = sd(m+1:M);
%         Tsph = mean(negv)/((prod(negv))^(1/(M-m)));
%         a(m+1) = 2*N*(M-m)*log(Tsph) + m*(2*M-m)*2;
%     end
%     [y,b] = min(a);
%     num = b - 1;
% else
%     [u,s,v] = svd(R);
%     sd = diag(s);
%     sd=sort(sd,'descend');
%     a = zeros(1,prod(M));
%     for m = 0:prod(M)-1
%         negv = sd(m+1:prod(M));
%         Tsph = mean(negv)/((prod(negv))^(1/(prod(M)-m)));
%         a(m+1) = N*(prod(M)-m)*log(Tsph) + m*(2*prod(M)-m);
%     end
%     [y,b] = min(a);
%     num = b - 1;
% end
% end