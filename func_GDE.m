%%
function  [GDE,Ns_GDE] = func_GDE(M,L,R)
Ns_GDE = 0;
Rp=R(1:M-1,1:M-1); 
% [U,D1]=eig(Rp); 
[U,S,V]=svd(Rp);
% U=fliplr(U); 
% T=zeros(M,M); 
% T(1:M-1,1:M-1)=U;
% T(M,M)=1; 
T=[U,zeros(M-1,1);zeros(1,M-1),1];
Rt=T'*R*T; 
r=abs(Rt(1:M-1,M)); 
% DL=exp(-L); 
DL = 0.3;
GDE=r-mean(r)*DL; 
for k=1:M-1
    if GDE(k)<0   
        Ns_GDE=k-1; break;
    end
end
% K=find(GDE>=0);   % �ҳ���һ���������������Ӧ��λ��
% Ns_GDE=M-K(1)+1;           % ��Դ��Ŀ�Ĺ���
% K=find(GDE>=0);%%%%%%%%%�ҳ���һ���������������Ӧ��λ��
% Ns_GDE=M-(K(1)-1);       %%%%%%%%%��Դ��Ŀ�Ĺ���
end

%%
% function [GDE,num]=func_GDE(N,R)
% %�Ƕ�Բ׼�������Դ��Ŀ
% %input R Э�������
% %input N  ������
% %output num ��Դ��
%     [row_length,col_length]=size(R);
%     R_new=R(1:row_length-1,1:col_length-1);
%     [V,D]=eig(R_new);
%     D =diag(D).';
%     [D,I0] = sort(D);
%     D=fliplr(D);
%     V=fliplr(V(:,I0));
%     U=[V zeros(row_length-1,1);zeros(1,col_length-1) 1];
%     S=U'*R*U;
%     k=1;
%     while 1
%         GDE=abs(S(k,col_length))-1/(2*N*(col_length-1))*sum(abs(S(1:row_length-1,col_length)));%��������ȡΪ1/N
%         if GDE<0 ||k>col_length-2
%             break;
%         end
%         k=k+1;
%     end
%     num=k-1;
% end