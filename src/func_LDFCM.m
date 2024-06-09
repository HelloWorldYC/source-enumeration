%% 论文复现
% 论文：A Source Number Estimation Algorithm Based on Data Local Density and Fuzzy C-Means Clustering
% 
%%
function [Ns_LDFCM]=func_LDFCM(R)

% L=100;
% num=3; %信源数
% Array_Num=10;% 阵元数
% d=0.5; %线阵半径
% lamda=1; %波长
% kk=6;    %线阵
% M=Array_Num;
% theta_jam1=10;
% theta_jam2=40;
% theta_jam3=70;
% alfa_jam1=10;
% alfa_jam2=50;
% alfa_jam3=90;
% s_jam1=array_form(Array_Num,d,lamda,theta_jam1,alfa_jam1,kk);
% s_jam2=array_form(Array_Num,d,lamda,theta_jam2,alfa_jam2,kk);
% s_jam3=array_form(Array_Num,d,lamda,theta_jam3,alfa_jam3,kk);
% A=[s_jam1;s_jam2;s_jam3];%方向矩阵；
% A=A';
% snr=-3;
% Am=10^(snr/10);
% x=randn(num,L);
% signal=Am*x;
% A1=A*signal; 
% X=awgn(A1,snr,'measured');
% R=X*X'/L;

[M,L]=size(R);
[U,S,V]=svd(R);
S = diag(S);
S = S+sqrt(sum(S));
s_min = min(S);
s_max = max(S);

dkl = zeros(M,M);
roukl = zeros(M,M);
rouk = zeros(1,M);
dc = sum(S)/M;

for k=1:M
    for l=1:M
        dkl(k,l)= sqrt((S(k)-S(l))^2);
        if (dkl(k,l)-dc) >= 0
            roukl(k,l) = 0;
        else
            roukl(k,l) = 1;
        end
    end
    roukl(k,k) = 0;
    rouk(k) = sum(roukl(k,:));
end

input = zeros(M,2);
for i=1:M
    input(i,1) = rouk(i);
    input(i,2) = S(i);
end

options = [2 200 1e-6 0];
[center, U, obj_fcn] = fcm(input, 2, options);

% figure;
% scatter(center(1,1),center(1,2),'rx');
% hold on;
% scatter(center(2,1),center(2,2),'bx');
% scatter(input(:,1),input(:,2),'go');

[max_U,index] = max(U);
Ns_LDFCM = length(find(index==index(1)));

end