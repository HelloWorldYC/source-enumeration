%%
function [MDL,Ns_MDL] = func_MDL(M,L,D)
%   MDL函数的实现
MDL=zeros(1,M-1);
for n=1:M-1
    temp0=sum(D(n+1:M))/(M-n);
    temp1=prod(D(n+1:M))^(1/(M-n));
    AIta=temp0/temp1;
    MDL(n)=L*(M-n)*log10(AIta)+(1/2)*n*(2*M-n)*log10(L);
end
[MDL_min,Ns_MDL] = min(MDL);
% Ns_MDL = M-(Ns_MDL-1)+1;
% Ns_MDL = Ns_MDL-1;
end

%%
% function [num]=func_MDL(R,N,M,D)
% %MDL准则估计信源数目
% %input R 阵元接收信号的协方差矩阵
% %input N 采样数
% %input M 阵元数
% %input D 1：一维线阵 2：二维方阵
% %output num 信源数目
% if D==1
%     [u,s,v] = svd(R);
%     sd = diag(s);
%     sd=sort(sd,'descend');%从大到小排列
%     a = zeros(1,M);
%     for m = 0:M-1
%         negv = sd(m+1:M);
%         Tsph = mean(negv)/((prod(negv))^(1/(M-m)));
%         a(m+1) = N*(M-m)*log(Tsph) + m*(2*M-m)*log(N)/2;
%     end
%     [y,b] = min(a);
%     num = b - 1;
% else
%     [u,s,v] = svd(R);
%     sd = diag(s);
%     sd=sort(sd,'descend');%从大到小排列
%     a = zeros(1,prod(M)); %prod([M,N])=M*N
%     for m = 0:prod(M)-1
%         negv = sd(m+1:prod(M));
%         Tsph = mean(negv)/((prod(negv))^(1/(prod(M)-m)));
%         a(m+1) = N*(prod(M)-m)*log(Tsph) + m*(2*prod(M)-m)*log(N)/2;
%     end
%     [y,b] = min(a);
%     num = b - 1;
% end
% end