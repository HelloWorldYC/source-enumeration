function [ DMDL ] = func_DMDL(M,L,D)
%   DMDL函数的实现
DMDL=zeros(1,M-1);
lamda_dl=sqrt(sum(D(1:M)));
mu=D+lamda_dl;
 for n=1:M-1
        temp0=sum(mu(n+1:M))/(M-n);
        temp1=1;
        for ii=n+1:M
            temp1=temp1*mu(ii);
        end
        temp1=temp1^(1/(M-n));
        AIta=temp0/temp1;
        %%%%%%%%%%%%%%%%%%
%         RAIC(n)=2*L*(M-n)*log(AIta)+2*n*(2*M-n);
          DMDL(n)=L*(M-n)*log(AIta)+1/2*n*(2*M-n)*log(L);
 end
end
