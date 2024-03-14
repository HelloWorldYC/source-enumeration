%% 程序说明
%  功能： 利用双树复小波分解单通道信号，再利用传统算法估计信源
%  码猿： 马叶椿
%  版本： v1.0 - 2021.12.16

%% 构造单通道信号
% 单通道信号构造比较简单，不需要通过方向矢量矩阵来构造，因为不存在阵元之间的关系
% 直接将几个信号源按增益叠加在一起即可，这个增益就是天线对信号源的增益，已经包括
% 了入射方向的影响
clear;
clc;
close all;
tic;

fs = 5000;
L = 1024;
num = 3;    % 信源数
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

s_jam1=array_form(1,d,lamda,theta_jam1,alfa_jam1,kk);
s_jam2=array_form(1,d,lamda,theta_jam2,alfa_jam2,kk);
s_jam3=array_form(1,d,lamda,theta_jam3,alfa_jam3,kk);
A=[s_jam1;s_jam2;s_jam3];%方向矩阵；
A=A';

%% 不同信噪比下的 Monte-Carlo 实验
monte=50;           %Monte-Carlo模拟的次数
biort = 'near_sym_b';
qshift = 'qshift_d';
nlevel = 10:1:10;
SNR = 10:1:30;
accuracy_SNR_AIC = zeros(length(nlevel),length(SNR));
accuracy_SNR_MDL = zeros(length(nlevel),length(SNR));
accuracy_SNR_GDE = zeros(length(nlevel),length(SNR));
accuracy_SNR_BIC = zeros(length(nlevel),length(SNR));
accuracy_SNR_AIC_EMD = zeros(1,length(SNR));
accuracy_SNR_MDL_EMD = zeros(1,length(SNR));
accuracy_SNR_GDE_EMD = zeros(1,length(SNR));
accuracy_SNR_BIC_EMD = zeros(1,length(SNR));
for level = 1:1:length(nlevel)
    unit = eye(nlevel(level));
    z = zeros(L,nlevel(level)+1);
    for i = 1:1:length(SNR)
        Am=10^(SNR(i)/10);
        Ns_AIC=zeros(2,monte);
        Ns_MDL=zeros(2,monte);
        Ns_GDE=zeros(2,monte);
        Ns_BIC=zeros(2,monte);
        for mk = 1:1:monte
            % 窄带高斯信号
%             [t1,at1,bt1,s1]=narrow_signal(fs,L,100,120,1500);
%             [t2,at2,bt2,s2]=narrow_signal(fs,L,80,95,2000);
%             [t3,at3,bt3,s3]=narrow_signal(fs,L,130,120,1800);
            % 低频随机信号
            x = randn(num,L);
            s1 = x(1,:);
            s2 = x(2,:);
            s3 = x(3,:);
            
            s = s1+s2+s3;
            noise = randn(1,L);
%             X = awgn(s,SNR(i),'measured');
            X = Am*s+noise;
            X = X';
           %% DCTWT分解
            %（二维的双树复小波有六个分解方向 ±15，±45，±75，但不意味着分解六层，
            % 对于一维来说，是没有这么多分解方向的，所以分解层数可选）
            [Yl,Yh,Yscale] = dtwavexfm(X,nlevel(level),biort,qshift);
            for j = 1:nlevel(level)
                z(:,j) = dtwaveifm(Yl*0,Yh,biort,qshift,unit(j,:));
            end
            z(:,nlevel(level)+1) = dtwaveifm(Yl,Yh,biort,qshift,zeros(1,nlevel(level)));
%             z = fliplr(z);
%             z(:,1:nlevel(level)) = z(:,1:nlevel(level))+z(:,nlevel(level)+1);
            Y1 = [z X];
            [L,M1] = size(Y1);
            Y1= Y1';
            R1 = Y1*Y1'/L;
            [u1,v1] = svd(R1);
            T1 = diag(v1);
            Tav1 = sqrt(sum(T1));
            T1 = T1+Tav1;
            [AIC,Ns_AIC(1,mk)] = func_AIC(M1,L,T1);
            [MDL,Ns_MDL(1,mk)] = func_MDL(M1,L,T1);
            [GDE,Ns_GDE(1,mk)] = func_GDE(M1,L,R1);
            [BIC,Ns_BIC(1,mk)] = func_NBIC(1/(M1*L),M1,L,R1);
           %% EMD分解
            emd_num = 8;
            Y2 = emd(X,'MaxNumIMF',emd_num);
            [M2,L]=size(Y2);
%             Y21 = Y2(1:M2-1,:)+Y2(M2,:);
            Y22 = [Y2 X];
            Y22 = Y22';
            [M2,L]=size(Y22);
            R2 = Y22*Y22'/L;
            [u2,v2] = svd(R2);
            T2 = diag(v2);
            Tav2 = sqrt(sum(T2));
            T2=T2+Tav2;
            [AIC,Ns_AIC(2,mk)] = func_AIC(M2,L,T2);
            [MDL,Ns_MDL(2,mk)] = func_MDL(M2,L,T2);
            [GDE,Ns_GDE(2,mk)] = func_GDE(M2,L,R2);
            [BIC,Ns_BIC(2,mk)] = func_IBIC(1/(M2*L),M2,L,R2);
        end
        accuracy_SNR_GDE(level,i)=length(find(Ns_GDE(1,:)==num))./monte;
        accuracy_SNR_MDL(level,i)=length(find(Ns_MDL(1,:)==num))./monte;
        accuracy_SNR_AIC(level,i)=length(find(Ns_AIC(1,:)==num))./monte;
        accuracy_SNR_BIC(level,i)=length(find(Ns_BIC(1,:)==num))./monte;
        accuracy_SNR_GDE_EMD(i)=length(find(Ns_GDE(2,:)==num))./monte;
        accuracy_SNR_MDL_EMD(i)=length(find(Ns_MDL(2,:)==num))./monte;
        accuracy_SNR_AIC_EMD(i)=length(find(Ns_AIC(2,:)==num))./monte;
        accuracy_SNR_BIC_EMD(i)=length(find(Ns_BIC(2,:)==num))./monte;
    end
end
%% 
figure(3);
plot(SNR,accuracy_SNR_GDE,'b^-');
hold on;
plot(SNR,accuracy_SNR_MDL,'gX-');
plot(SNR,accuracy_SNR_AIC,'o-');
plot(SNR,accuracy_SNR_BIC,'p-');
plot(SNR,accuracy_SNR_GDE_EMD,'ks-');
plot(SNR,accuracy_SNR_MDL_EMD,'r*-');
plot(SNR,accuracy_SNR_AIC_EMD,'c*-');
plot(SNR,accuracy_SNR_BIC_EMD,'g*-');
xlabel('SNR');
ylabel('accuracy');
title(['DTCWT分解层数为:',num2str(nlevel),'    EMD分解层数为:',num2str(emd_num)]);
legend({'DTCWT_GDE','DTCWT_MDL','DTCWT_AIC','DTCWT_BIC','EMD_GDE','EMD_MDL','EMD_AIC','EMD_BIC'},'Interpreter','none');

toc;