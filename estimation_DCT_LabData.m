%% ����˵��
%  ���ܣ� ����˫����С���ֽⵥͨ���źţ�ʵ���ҵ����ݣ��������ô�ͳ�㷨������Դ
%  ��Գ�� ��Ҷ��
%  �汾�� v1.0 - 2022.01.14

%% ���뵥ͨ���ź�

clear;
clc;
close all;
tic;

name1 = 'EW2_190_z12_waveform';
num = 2;
snr = 12;
x1 = xlsread(['F:\��Դ������\code\MYC\DTCWT\labdata\',name1],'D1100:D54000');
L_total = length(x1);

%% ��ͬ������µ� Monte-Carlo ʵ��
monte=50;           %Monte-Carloģ��Ĵ���
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
    [Yl,Yh,Yscale] = dtwavexfm(signal,nlevel,biort,qshift);%Yh�дӸ�Ƶ����Ƶ
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
title(['��Դ����',num2str(num),'    ����ȣ�',num2str(snr),'    DCTWT�ֽ������',num2str(nlevel)]);
legend(['GDE��',num2str(accuracy_GDE)]);

toc;