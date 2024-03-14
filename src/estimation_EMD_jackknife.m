%% 程序说明
%  功能： 先利用EMD分解单通道信号，再用刀切法切段，再来估计
%  码猿： 马叶椿
%  版本： v1.0 - 2022.01.16
%         v2.0 - 2022.01.18
%         v3.0 - 2022.05.14
%% 输入单通道信号

clear;
clc;
close all;
tic;

f0 = 15e6;
fs = 60e6;
fa = 10e6;
fb = 10e6;
L = 512;
num = 3;

Array_Num=1;% 阵元数
d=0.5; %线阵半径
lamda=1; %波长
kk=6;    %线阵
% 入射角
theta_jam1=20;
theta_jam2=30;
theta_jam3=45;
degrad=pi/180;
%方位角
alfa_jam1=10;
alfa_jam2=50;
alfa_jam3=90;

s_jam1=array_form(Array_Num,d,lamda,theta_jam1,alfa_jam1,kk);
s_jam2=array_form(Array_Num,d,lamda,theta_jam2,alfa_jam2,kk);
s_jam3=array_form(Array_Num,d,lamda,theta_jam3,alfa_jam3,kk);
A=[s_jam1;s_jam2;s_jam3];%方向矩阵；
A=A';

Nt=100; %Monte次数
jj=0;
rou = 0.75;
times = 40;
snr=-20:20;
Pd_GDE=zeros(1,length(snr));
Pd_AIC=zeros(1,length(snr));
Pd_MDL=zeros(1,length(snr));
Pd_BIC=zeros(1,length(snr));
Pd_MIC=zeros(1,length(snr));
Pd_MSTDC=zeros(1,length(snr));


for SNR=snr 
    disp(['SNR is ',num2str(SNR)]);
    source_power=10.^(SNR./10);
    source_amplitude = sqrt(source_power)*ones(1,num);    % 信源标准差
    jj=jj+1;
    Ns_AIC=zeros(1,Nt);
    Ns_MDL=zeros(1,Nt);
    Ns_GDE=zeros(1,Nt);
    Ns_BIC=zeros(1,Nt);
    Ns_MIC=zeros(1,Nt);
    Ns_MSTDC=zeros(1,Nt);

for cc=1:Nt
    frequency_GDE=zeros(1,times);
    frequency_AIC=zeros(1,times);
    frequency_MDL=zeros(1,times);
    frequency_BIC=zeros(1,times);
    frequency_MIC=zeros(1,times);
    frequency_MSTDC=zeros(1,times);

    [t1,at1,bt1,x1]=narrow_signal(fs,L,fa,fb,f0);
    [t2,at2,bt2,x2]=narrow_signal(fs,L,fa,fb,f0);
    [t3,at3,bt3,x3]=narrow_signal(fs,L,fa,fb,f0);
%     [code1,x1] = BPSK_generate(500,1000,2000,128,1);
%     [code2,x2] = BPSK_generate(500,1000,2000,128,1);
%     [code3,x3] = BPSK_generate(500,1000,2000,128,1);
    x = [x1;x2;x3];
    source_wave = sqrt(0.5)*x;
    st=source_amplitude*source_wave;
    signal = awgn(st,SNR,'measured');
    signal = signal';
    
    z = emd(signal);
    Y = [signal z];
    Y = Y';
    [M,~] = size(Y);
    [signal_jk,index,L1] = Jackknife(Y(1,:),rou,times);
    subsets = cell(times,1);
    subsets(:) = {zeros(M,L1)};
    for i=1:times
        for j=1:M
            subsets{i}(j,:) = Y(j,index(i,:));
        end
    end
    for k=1:times
        resignal = subsets{k};
        R = resignal*resignal'/L1;
        [~,v]=svd(R);
        T=diag(v);
        T1=T+sqrt(sum(T));
        [AIC,frequency_AIC(k)] = func_AIC(M,L1,T1);
        [MDL,frequency_MDL(k)] = func_MDL(M,L1,T1);
        [GDE,frequency_GDE(k)] = func_GDE(M,L1,R);
        [BIC,frequency_BIC(k)] = func_NBIC(1/(M*L1),M,L1,R);
        [MIC,frequency_MIC(k)] = func_MIC(resignal,M,L1);
        [MSTDC,frequency_MSTDC(k)] = func_MSTDC(resignal,M,L1);
    end
    Ns_GDE(cc) = mode(frequency_GDE);
    Ns_AIC(cc) = mode(frequency_AIC);
    Ns_MDL(cc) = mode(frequency_MDL);
    Ns_BIC(cc) = mode(frequency_BIC);
    Ns_MIC(cc) = mode(frequency_MIC);
    Ns_MSTDC(cc) = mode(frequency_MSTDC);
  
end

Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt;
Pd_MDL(jj)=length(find(Ns_MDL==num))./Nt;
Pd_AIC(jj)=length(find(Ns_AIC==num))./Nt;
Pd_BIC(jj)=length(find(Ns_BIC==num))./Nt;
Pd_MIC(jj)=length(find(Ns_MIC==num))./Nt;
Pd_MSTDC(jj)=length(find(Ns_MSTDC==num))./Nt;

end
%%
plot(snr,Pd_GDE,'>-',snr,Pd_BIC,'rs-',snr,Pd_MIC,'o-',snr,Pd_MSTDC,'b^-',snr,Pd_MDL,'gs-');
% plot(xx,Pd_AIC,'>-',xx,Pd_MDL,'rs-',xx,Pd_GDE,'o-',xx,Pd_BIC,'b^-',xx,Pd_MIC,'c*-',xx,Pd_MSTDC,'r*-');
title(['emd分解估计',num2str(num),'个信源']);
xlabel('不同信噪比（dB）');
ylabel('正确检测概率(%)');
axis([min(snr) max(snr) 0 1]);
legend('EJGDE','EJBIC','EJAMIC','EJAMSTDC','EJDMDL');
toc;