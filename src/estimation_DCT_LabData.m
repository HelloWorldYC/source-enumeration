%% 程序说明
%  功能： 利用双树复小波分解单通道信号（实验室的数据），再利用传统算法估计信源
%  码猿： 马叶椿
%  版本： v1.0 - 2022.01.14

%% 输入单通道信号

clear;
clc;
close all;
tic;

name1 = 'EW2_190_z12_waveform';
num = 2;
snr = 12;
x1 = xlsread(['F:\信源数估计\code\MYC\DTCWT\labdata\',name1],'D1100:D54000');
L_total = length(x1);

%% 不同信噪比下的 Monte-Carlo 实验
monte=50;           %Monte-Carlo模拟的次数
L = 1024;
% biort = 'antonini';
biort = 'near_sym_b';
qshift = 'qshift_b';
nlevel = 4;
unit = eye(nlevel);
z = zeros(L,nlevel+1);
Ns_GDE = zeros(1,monte);
for mk = 1:1:monte
    index = (mk-1)*L+1:mk*L;
    signal = x1(index);
    [Yl,Yh,Yscale] = dtwavexfm(signal,nlevel,biort,qshift);%Yh中从高频到低频
    for j = 1:nlevel
        z(:,j) = dtwaveifm(Yl,Yh,biort,qshift,unit(j,:));
    end
%     z(:,nlevel+1) = dtwaveifm(Yl,Yh,biort,qshift,zeros(1,nlevel));
%     z(:,1:nlevel)=z(:,1:nlevel)+z(:,nlevel+1);
%         z = fliplr(z);
    Y = [z(:,1:nlevel) signal];
    Y = Y';
    [M,L] = size(Y);
    R = Y*Y'/L;
    [GDE,Ns_GDE(mk)]=func_GDE(M,L,R);

end
accuracy_GDE=length(find(Ns_GDE==num))./monte;

%% 
figure(3);
plot((1:monte),Ns_GDE,'b^-');
hold on;
% plot(SNR,accuracy_SNR_AIC(1,:),'gX-');
% plot(SNR,accuracy_SNR_MDL(1,:),'o-');
xlabel('SNR');
ylabel('estimation');
title(['信源数：',num2str(num),'    信噪比：',num2str(snr),'    DCTWT分解层数：',num2str(nlevel)]);
legend(['GDE：',num2str(accuracy_GDE)]);

toc;