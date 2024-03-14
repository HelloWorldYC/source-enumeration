%% 程序说明
%  功能�? 利用双树复小波分解单通道信号，再利用传统算法估计信源
%  码猿�? 马叶�?
%  版本�? v1.0 - 2021.12.16

%% 构�?�单通道信号
% 单�?�道信号构�?�比较简单，不需要�?�过方向矢量矩阵来构造，因为不存在阵元之间的关系
% 直接将几个信号源按增益叠加在�?起即可，这个增益就是天线对信号源的增益，已经包括
% 了入射方向的影响
clear;
clc;
close all;
tic;

f0 = 15e6;
fs = 62e6;
fa = 8.5e6;
fb = 8.5e6;
Ns = 512;
L=Ns;
num = 3;    % 信源�?
d=0.5; %线阵半径
lamda=1; %波长
kk=6;    %线阵

% 入射�?
theta_jam1=10;
theta_jam2=40;
theta_jam3=70;

%方位�?
alfa_jam1=10;
alfa_jam2=50;
alfa_jam3=90;

s_jam1=array_form(1,d,lamda,theta_jam1,alfa_jam1,kk);
s_jam2=array_form(1,d,lamda,theta_jam2,alfa_jam2,kk);
s_jam3=array_form(1,d,lamda,theta_jam3,alfa_jam3,kk);
A=[s_jam1;s_jam2;s_jam3];%方向矩阵�?
A=A';

%% 不同信噪比下�? Monte-Carlo 实验
monte=50;           %Monte-Carlo模拟的次�?
biort = 'near_sym_b';
qshift = 'qshift_d';
nlevel = 10;
SNR = -5:1:20;
rou = 0.75;
times = 20;
accuracy_SNR_AIC = zeros(1,length(SNR));
accuracy_SNR_MDL = zeros(1,length(SNR));
accuracy_SNR_GDE = zeros(1,length(SNR));
accuracy_SNR_NBIC = zeros(1,length(SNR));
accuracy_SNR_AIC_EMD = zeros(1,length(SNR));
accuracy_SNR_MDL_EMD = zeros(1,length(SNR));
accuracy_SNR_GDE_EMD = zeros(1,length(SNR));
accuracy_SNR_BIC_EMD = zeros(1,length(SNR));
accuracy_SNR_jackkingfe_emd_MIC = zeros(1,length(SNR));
accuracy_SNR_jackkingfe_emd_MSTDC = zeros(1,length(SNR));

unit = eye(nlevel);

parfor i = 1:1:length(SNR)
    disp(['SNR is ',num2str(SNR(i))]);
    L = 512;
    z = zeros(L,nlevel+1);
    Am=10^(SNR(i)/10);
    Ns_AIC=zeros(2,monte);
    Ns_MDL=zeros(2,monte);
    Ns_GDE=zeros(2,monte);
    Ns_NBIC=zeros(2,monte);
    Ns_jackknife_emd_MIC=zeros(1,monte);
    Ns_jackknife_emd_MSTDC=zeros(1,monte);
    for mk = 1:1:monte
        % 窄带高斯信号
%       [t1,at1,bt1,s1]=narrow_signal(fs,L,100,120,1500);
%       [t2,at2,bt2,s2]=narrow_signal(fs,L,80,95,2000);
%       [t3,at3,bt3,s3]=narrow_signal(fs,L,130,120,1800);
        [t1,at1,bt1,s1]=narrow_signal(fs,L,fa,fb,f0);
        [t2,at2,bt2,s2]=narrow_signal(fs,L,fa,fb,f0);
        [t3,at3,bt3,s3]=narrow_signal(fs,L,fa,fb,f0);
        % 低频随机信号
%         x = randn(num,L);
%         s1 = x(1,:);
%         s2 = x(2,:);
%         s3 = x(3,:);
            
        s = 0.5*s1+0.5*s2+0.5*s3;
        noise = randn(1,L);
%         X = awgn(s,SNR(i),'measured');
        X = Am*s+noise;
        X = X';
        [Ns_jackknife_emd_MIC(mk),Ns_jackknife_emd_MSTDC(mk)] = func_jackkinfe_emd(X,rou,times);
       %% DCTWT分解
        %（二维的双树复小波有六个分解方向 ±15，�?45，�?75，但不意味着分解六层�?
        % 对于�?维来说，是没有这么多分解方向的，�?以分解层数可选）
        [Yl,Yh,Yscale] = dtwavexfm(X,nlevel,biort,qshift);
        
%         for j = 1:nlevel
%             z(:,j) = dtwaveifm(Yl*0,Yh,biort,qshift,unit(j,:));
%         end
        z(:,1) = dtwaveifm(Yl*0,Yh,biort,qshift,unit(1,:));
        z(:,2) = dtwaveifm(Yl*0,Yh,biort,qshift,unit(2,:));
        z(:,3) = dtwaveifm(Yl*0,Yh,biort,qshift,unit(3,:));
        z(:,4) = dtwaveifm(Yl*0,Yh,biort,qshift,unit(4,:));
        z(:,5) = dtwaveifm(Yl*0,Yh,biort,qshift,unit(5,:));
        z(:,6) = dtwaveifm(Yl*0,Yh,biort,qshift,unit(6,:));
        z(:,7) = dtwaveifm(Yl*0,Yh,biort,qshift,unit(7,:));
        z(:,8) = dtwaveifm(Yl*0,Yh,biort,qshift,unit(8,:));
        z(:,9) = dtwaveifm(Yl*0,Yh,biort,qshift,unit(9,:));
        z(:,10) = dtwaveifm(Yl*0,Yh,biort,qshift,unit(10,:));

        z(:,nlevel+1) = dtwaveifm(Yl,Yh,biort,qshift,zeros(1,nlevel));
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
        [NBIC,Ns_NBIC(1,mk)] = func_NBIC(1/(M1*L),M1,L,R1);
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
         [BIC,Ns_NBIC(2,mk)] = func_IBIC(1/(M2*L),M2,L,R2);
    end
    accuracy_SNR_GDE(i)=length(find(Ns_GDE(1,:)==num))./monte;
    accuracy_SNR_MDL(i)=length(find(Ns_MDL(1,:)==num))./monte;
    accuracy_SNR_AIC(i)=length(find(Ns_AIC(1,:)==num))./monte;
    accuracy_SNR_NBIC(i)=length(find(Ns_NBIC(1,:)==num))./monte;
    accuracy_SNR_GDE_EMD(i)=length(find(Ns_GDE(2,:)==num))./monte;
    accuracy_SNR_MDL_EMD(i)=length(find(Ns_MDL(2,:)==num))./monte;
    accuracy_SNR_AIC_EMD(i)=length(find(Ns_AIC(2,:)==num))./monte;
    accuracy_SNR_BIC_EMD(i)=length(find(Ns_NBIC(2,:)==num))./monte;
    accuracy_SNR_jackkingfe_emd_MIC(i)=length(find(Ns_jackknife_emd_MIC==num))./monte;
    accuracy_SNR_jackkingfe_emd_MSTDC(i)=length(find(Ns_jackknife_emd_MSTDC==num))./monte;
end
%% 
figure(3);
% plot(SNR,accuracy_SNR_GDE,'b^-');
hold on;
plot(SNR,accuracy_SNR_MDL,'gX-');
plot(SNR,accuracy_SNR_AIC,'o-');
plot(SNR,accuracy_SNR_NBIC,'p-');
plot(SNR,accuracy_SNR_GDE_EMD,'ks-');
plot(SNR,accuracy_SNR_MDL_EMD,'r*-');
plot(SNR,accuracy_SNR_AIC_EMD,'c*-');
plot(SNR,accuracy_SNR_BIC_EMD,'g*-');
plot(SNR,accuracy_SNR_jackkingfe_emd_MIC,'rs-');
plot(SNR,accuracy_SNR_jackkingfe_emd_MSTDC,'gs-');
xlabel('SNR');
ylabel('accuracy');
% title(['DTCWT分解层数�?:',num2str(nlevel),'    EMD分解层数�?:',num2str(emd_num)]);
legend({'DTCWT_MDL','DTCWT_AIC','DTCWT_BIC','EMD_GDE','EMD_MDL','EMD_AIC','EMD_BIC','jk_emd_MIC','jk_emd_MSTDC'},'Interpreter','none');

%%
% f0 = 15e6;
% fs = 62e6;
% fa = 3e6;
% fb = 3e6;
% Ns = 1024;
% L=Ns;
[t1,at1,bt1,s1]=narrow_signal(fs,L,fa,fb,f0);
[t2,at2,bt2,s2]=narrow_signal(fs,L,fa,fb,f0);
[t3,at3,bt3,s3]=narrow_signal(fs,L,fa,fb,f0);
s = s1+s2;
noise = randn(1,L);
Am=10^(15/10);
X = Am*s+noise;
% X = awgn(s,20,'measured');
Rtau=xcorr(X);                            %自相关函�?
Sx=fft(Rtau);                             %功率谱密�?
len=length(Sx);
k=0:len-1;
w=2*pi*(k/len-1/2)*fs;
figure,plot(w/2/pi,abs(fftshift(Sx)));
title('功率谱密度S_x(\omega)'),xlabel('f/Hz');

toc;