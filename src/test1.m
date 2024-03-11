%% 生成阵列多通道信号，检测AIC、MDL、GDE、BIC的性能
clear;
clc;

tic;
f0 = 15.48e4;
fs = 62e4;
fa = 2.3e4;
fb = 2.2e4;
Ns = 512;
L=Ns;
t=1:L;
num=2; %信源数

Array_Num=10;% 阵元数
d=0.5; %线阵半径
lamda=1; %波长
kk=6;    %线阵
M=Array_Num;
% 入射角
theta_jam1=10;
theta_jam2=40;
theta_jam3=70;
theta_jam4=70;
theta_jam5=88;
degrad=pi/180;
%方位角
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
A=[s_jam1;s_jam2];%方向矩阵；
% A=[s_jam1];
A=A';

Nt=100; %Monte次数
jj=0;
Pd_GDE=zeros(1,41);
Pd_AIC=zeros(1,41);
Pd_MDL=zeros(1,41);
Pd_IBIC=zeros(1,41);
Pd_MIC=zeros(1,41);
Pd_MSTDC=zeros(1,41);
Pd_ISSM=zeros(1,41);
Pd_LDFCM=zeros(1,41);
for SNR=-20:20 
    disp(['SNR is ',num2str(SNR)]);
    Am=10^(SNR/10);
    jj=jj+1;
    Ns_AIC=zeros(1,Nt);
    Ns_MDL=zeros(1,Nt);
    Ns_GDE=zeros(1,Nt);
    Ns_IBIC=zeros(1,Nt);
    Ns_MIC=zeros(1,Nt);
    Ns_MSTDC=zeros(1,Nt);
    Ns_ISSM=zeros(1,Nt);
    Ns_LDFCM=zeros(1,Nt);
for cc=1:Nt
    [t1,at1,bt1,x1]=narrow_signal(fs,L,fa,fb,f0);
    [t2,at2,bt2,x2]=narrow_signal(fs,L,fa,fb,f0);
    [t3,at3,bt3,x3]=narrow_signal(fs,L,fa,fb,f0);
    x = [x1;x2];
%     x=randn(num,L);
    signal=Am*x;
    A1=A*signal; 
    X=awgn(A1,SNR,'measured');
    X = real(X);
%     noise=randn(M,L); %白噪声模型
%     % noise=randn(M,L); %白噪声模型
%     X=A1+noise;

    R=X*X'/L; %信号协方差
    [u,v]=svd(R);
    T=diag(v);
    % T1=T+sqrt(sum(T));
    [AIC,Ns_AIC(cc)] = func_AIC(M,L,T);
    [MDL,Ns_MDL(cc)] = func_MDL(M,L,T);
    [GDE,Ns_GDE(cc)] = func_GDE(M,L,R);
    [IBIC,Ns_IBIC(cc)] = func_IBIC(1/(M*L),M,L,R);
    [MIC,Ns_MIC(cc)] = func_MIC(X,M,L);
    [alpha,Ns_MSTDC(cc)] = func_MSTDC(X,M,L);
    [ISSM,Ns_ISSM(cc)]=func_ISSM(X);
    [Ns_LDFCM(cc)]=func_LDFCM(R);
end

Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt;
Pd_MDL(jj)=length(find(Ns_MDL==num))./Nt;
Pd_AIC(jj)=length(find(Ns_AIC==num))./Nt;
Pd_IBIC(jj)=length(find(Ns_IBIC==num))./Nt;
Pd_MIC(jj)=length(find(Ns_MIC==num))./Nt;
Pd_MSTDC(jj)=length(find(Ns_MSTDC==num))./Nt;
Pd_ISSM(jj)=length(find(Ns_ISSM==num))./Nt;
Pd_LDFCM(jj)=length(find(Ns_LDFCM==num))./Nt;
end
 %%
xx=-20:20;
plot(xx,Pd_MDL,'rs-',xx,Pd_GDE,'o-',xx,Pd_IBIC,'b^-',xx,Pd_MIC,'gs-',xx,Pd_ISSM,'bs-',xx,Pd_LDFCM,'>-');
title(['白噪声下',num2str(Array_Num),'线阵估计',num2str(num),'个信源']);
xlabel('不同信噪比（dB）');
ylabel('正确检测概率(%)');
axis([-20 20 0 1]);
legend('MDL','GDE','IBIC','MIC','ISSM','LDFCM');
toc;