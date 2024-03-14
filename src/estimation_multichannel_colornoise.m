%% 多�?�道色噪声下估计
clear;
clc;
tic;

f0 = 15.48e4;
fs = 62e4;
fa = 2.3e4;
fb = 2.2e4;
Ns = 256;
L=Ns;
t=1:L;
num=3; %信源�?

Array_Num=8;% 阵元�?
d=0.5; %线阵半径
lamda=1; %波长
kk=6;    %线阵
M=Array_Num;
% 入射�?
theta_jam1=10;
theta_jam2=40;
theta_jam3=70;
theta_jam4=70;
theta_jam5=88;
degrad=pi/180;
%方位�?
alfa_jam1=10;
alfa_jam2=50;
alfa_jam3=90;
alfa_jam4=140;
alfa_jam5=190;

s_jam1=array_form(Array_Num,d,lamda,theta_jam1,alfa_jam1,kk);
s_jam2=array_form(Array_Num,d,lamda,theta_jam2,alfa_jam2,kk);
s_jam3=array_form(Array_Num,d,lamda,theta_jam3,alfa_jam3,kk);
s_jam4=array_form(Array_Num,d,lamda,theta_jam4,alfa_jam4,kk);
s_jam5=array_form(Array_Num,d,lamda,theta_jam5,alfa_jam5,kk);
A=[s_jam1;s_jam2;s_jam3];%方向矩阵�?
% A=[s_jam1];
A=A';

% 构�?�低通滤波器
Wp=2*pi*30;
Ws=2*pi*40;
Rp=0.5;
Rs=40;
fs1=100;
W=2*pi*fs1;
[N1,Wn]=buttord(2*Wp/W,2*Ws/W,Rp,Rs);
[b,a]=butter(N1,Wn);

Nt=50; %Monte次数
jj=0;
snr_min = -20;
snr_max = 20;
snr_length = snr_max-snr_min+1;
Pd_GDE=zeros(1,snr_length);
Pd_AIC=zeros(1,snr_length);
Pd_MDL=zeros(1,snr_length);
Pd_NBIC=zeros(1,snr_length);
Pd_MIC=zeros(1,snr_length);
Pd_MSTDC=zeros(1,snr_length);
Pd_ISSM=zeros(1,snr_length);
Pd_LDFCM=zeros(1,snr_length);

for SNR=snr_min:snr_max 
    disp(['SNR is ',num2str(SNR)]);
    Am=10^(SNR/10);
    jj=jj+1;
    Ns_AIC=zeros(1,Nt);
    Ns_MDL=zeros(1,Nt);
    Ns_GDE=zeros(1,Nt);
    Ns_NBIC=zeros(1,Nt);
    Ns_MIC=zeros(1,Nt);
    Ns_MSTDC=zeros(1,Nt);
    Ns_ISSM=zeros(1,Nt);
    Ns_LDFCM=zeros(1,Nt);

for cc=1:Nt
%     x=randn(num,L);
    [t1,at1,bt1,x1]=narrow_signal(fs,L,fa,fb,f0);
    [t2,at2,bt2,x2]=narrow_signal(fs,L,fa,fb,f0);
    [t3,at3,bt3,x3]=narrow_signal(fs,L,fa,fb,f0);
%     [code1,x1] = BPSK_generate(500,1000,2000,128,1);
%     [code2,x2] = BPSK_generate(500,1000,2000,128,1);
%     [code3,x3] = BPSK_generate(500,1000,2000,128,1);
    x = [x1;x2;x3];
    signal=Am*x;
    A1=A*signal; 
    noise=randn(M,L);
    color_noise=filter(b,a,noise);        %滤波产生高斯色噪�?
    X=A1+color_noise;
%     X=awgn(A1,SNR,'measured');
%     noise=randn(M,L); %白噪声模�?
%     % noise=randn(M,L); %白噪声模�?
%     X=A1+noise;

    R=X*X'/L; %信号协方�?

    [u,v]=svd(R);
    T=diag(v);
    % T1=T+sqrt(sum(T));
    [AIC,Ns_AIC(cc)] = func_AIC(M,L,T);
    [MDL,Ns_MDL(cc)] = func_MDL(M,L,T);
    [GDE,Ns_GDE(cc)] = func_GDE(M,L,R);
    [IBIC,Ns_NBIC(cc)] = func_NBIC(1/(M*L),M,L,R);
    [MIC,Ns_MIC(cc)] = func_MIC(X,M,L);
    [MSTDC,Ns_MSTDC(cc)] = func_MSTDC(X,M,L);
    [ISSM,Ns_ISSM(cc)]=func_ISSM(X);
    [Ns_LDFCM(cc)]=func_LDFCM(R);

end

Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt;
Pd_MDL(jj)=length(find(Ns_MDL==num))./Nt;
Pd_AIC(jj)=length(find(Ns_AIC==num))./Nt;
Pd_NBIC(jj)=length(find(Ns_NBIC==num))./Nt;
Pd_MIC(jj)=length(find(Ns_MIC==num))./Nt;
Pd_MSTDC(jj)=length(find(Ns_MSTDC==num))./Nt;
Pd_ISSM(jj)=length(find(Ns_ISSM==num))./Nt;
Pd_LDFCM(jj)=length(find(Ns_LDFCM==num))./Nt;

end
 %%
xx=snr_min:snr_max;
plot(xx,Pd_GDE,'o-',xx,Pd_NBIC,'b^-',xx,Pd_MIC,'gs-',xx,Pd_ISSM,'>-',xx,Pd_LDFCM,'bs-');
% plot(xx,Pd_AIC,'>-',xx,Pd_MDL,'rs-',xx,Pd_GDE,'o-',xx,Pd_BIC,'b^-',xx,Pd_MIC,'c*-',xx,Pd_MSTDC,'r*-');
title(['色噪声下',num2str(Array_Num),'线阵估计',num2str(num),'个信�?']);
xlabel('不同信噪比（dB�?');
ylabel('正确�?测概�?(%)');
axis([snr_min snr_max 0 1]);
legend('GDE','NBIC','MIC','ISSM','LDFCM');
toc;
