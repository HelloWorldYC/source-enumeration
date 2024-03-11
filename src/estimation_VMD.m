%% 程序说明
%  功能： 利用VMD分解单通道信号，再利用传统算法估计信源
%  码猿： 马叶椿
%  版本： v1.0 - 2022.01.17

%% 输入单通道信号

clear;
clc;
close all;
tic;

fs = 1024;
ts = 1/fs;
L = 1024;
t = (0:L-1)*ts;
num = 3;
d=0.5; %线阵半径
lamda=1; %波长
kk=6;    %线阵
% 入射角
theta_jam1=10;
theta_jam2=40;
theta_jam3=70;
%方位角
alfa_jam1=10;
alfa_jam2=50;
alfa_jam3=90;

s_jam1=array_form(8,d,lamda,theta_jam1,alfa_jam1,kk);
s_jam2=array_form(8,d,lamda,theta_jam2,alfa_jam2,kk);
s_jam3=array_form(8,d,lamda,theta_jam3,alfa_jam3,kk);
A=[s_jam1;s_jam2;s_jam3];%方向矩阵；
A=A';

%%
SNR = 20:20;
monte = 100;
jj = 0;
Pd_GDE=zeros(1,length(SNR));
Pd_AIC=zeros(1,length(SNR));
Pd_MDL=zeros(1,length(SNR));
Pd_BIC=zeros(1,length(SNR));
Ns_AIC = zeros(1,monte);
Ns_MDL = zeros(1,monte);
Ns_GDE = zeros(1,monte);
Ns_BIC = zeros(1,monte);

for snr = SNR
    jj = jj+1;
    Am=10^(snr/10);
    for i = 1:monte
        x = randn(num,L);
        signal=Am*x;
        A1=A*signal; 
        A2=A1(3,:);
        X = awgn(A2,snr,'measured');
        num_mode = 10;
        alpha = 1200;
        tau = 0;
        [Y, u_hat, omega] = VMD(X, alpha, tau, num_mode, false, 1, 1e-7);
        %%
        % Y = fliplr(Y);
%         [M,L]=size(Y);
%         Y1=zeros(num_mode/2,L);
%         for j=1:num_mode/2
%             Y1(j,:) = Y(j,:)+Y(num_mode-j+1,:);
%         end
        Y2 = Y;
        [M,L]=size(Y2);
        % 自相关系数矩阵
%         y_mean = mean(Y2,2);
%         for v=1:M
%            Y2(v,:)=Y2(v,:)-y_mean(v); 
%         end
%         V=Y2*Y2'/L;
%         V_diag=diag(diag(V));
%         R = (V_diag^(-1/2))*V*(V_diag^(-1/2));

        R = Y2*Y2'/L;
        [GDE,Ns_GDE(i)] = func_GDE(M,L,R);
        [BIC,Ns_BIC(i)] = func_BIC(1/(M*L),M,L,R);
        [u,v]=svd(R);
        T=diag(v);
        T1=T+sqrt(sum(T));
        [AIC,Ns_AIC(i)] = func_AIC(M,L,T1);
        [MDL,Ns_MDL(i)] = func_MDL(M,L,T1);
    end
%%
    Pd_GDE(jj)=length(find(Ns_GDE==num))./monte;
    Pd_MDL(jj)=length(find(Ns_MDL==num))./monte;
    Pd_AIC(jj)=length(find(Ns_AIC==num))./monte;
    Pd_BIC(jj)=length(find(Ns_BIC==num))./monte;
end
%%
figure(1);
plot(SNR,Pd_GDE,'b^-');
hold on;
plot(SNR,Pd_AIC,'gX-');
plot(SNR,Pd_MDL,'o-');
plot(SNR,Pd_BIC,'rs-');
xlabel('SNR');
ylabel('accuracy');
legend('GDE','AIC','MDL','BIC');
title(['信源个数为',num2str(num),'       分解层数为',num2str(num_mode)]);
toc;