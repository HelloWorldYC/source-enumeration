%% 多通道色噪声环境下，通过白化滤波器将其转为白噪声，再估计信源

clear;
clc;

Ns = 1024;
L=Ns;
t=1:L;
num=3; %信源数

Array_Num=8;% 阵元数
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
A=[s_jam1;s_jam2;s_jam3];%方向矩阵；
% A=[s_jam1];
A=A';

% 构造低通滤波器
Wp=2*pi*30;
Ws=2*pi*40;
Rp=0.5;
Rs=40;
fs=100;
W=2*pi*fs;
[N1,Wn]=buttord(2*Wp/W,2*Ws/W,Rp,Rs);
[b,a]=butter(N1,Wn);

Nt=100; %Monte次数
jj=0;
Pd_GDE=zeros(1,41);
Pd_AIC=zeros(1,41);
Pd_MDL=zeros(1,41);
Pd_BIC=zeros(1,41);
Pd_BIC_qita10=zeros(1,41);
Pd_BIC_qita100=zeros(1,41);
for SNR=-30:10 
    Am=10^(SNR/10);
    jj=jj+1;
    Ns_AIC=zeros(1,Nt);
    Ns_MDL=zeros(1,Nt);
    Ns_GDE=zeros(1,Nt);
    Ns_BIC=zeros(1,Nt);
    Ns_BIC_qita10=zeros(1,Nt);
    Ns_BIC_qita100=zeros(1,Nt);
for cc=1:Nt
    x=randn(num,L);
    signal=Am*x;
    A1=A*signal; 
    noise=randn(M,L);
    color_noise=filter(b,a,noise);        %滤波产生高斯色噪声
    X=A1+color_noise;

    R=X*X'/L; %源信号协方差
    %%%%%%%%%%%%%改动：构造白化滤波器%%%%%%%%%%%%%%%
    Sn=color_noise*color_noise'/L;%纯噪声协方差矩阵
    W=Sn^(-1/2);%构造白化滤波器
    Y=W*R*W';%白化滤波变换后的协方差矩阵
    %         Y=W*Rs*W';%白化滤波变换后的协方差矩阵
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
[u,v]=svd(Y);
T=diag(v);
% T1=T+sqrt(sum(T));
[AIC,Ns_AIC(cc)] = func_AIC(M,L,T);
[MDL,Ns_MDL(cc)] = func_MDL(M,L,T);
[GDE,Ns_GDE(cc)] = func_GDE(M,L,Y);
[BIC,Ns_BIC(cc)] = func_BIC(1/(M*L),M,L,Y);
[BIC_qita10,Ns_BIC_qita10(cc)] = func_BIC(10/(M*L),M,L,Y);
[BIC_qita100,Ns_BIC_qita100(cc)] = func_BIC(100/(M*L),M,L,R);
end

Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt;
Pd_MDL(jj)=length(find(Ns_MDL==num))./Nt;
Pd_AIC(jj)=length(find(Ns_AIC==num))./Nt;
Pd_BIC(jj)=length(find(Ns_BIC==num))./Nt;
Pd_BIC_qita10(jj)=length(find(Ns_BIC_qita10==num))./Nt;
Pd_BIC_qita100(jj)=length(find(Ns_BIC_qita100==num))./Nt;
end
 %%
figure(1);
xx=-30:10;
plot(xx,Pd_AIC,'>-',xx,Pd_MDL,'rs-',xx,Pd_GDE,'o-',xx,Pd_BIC,'b^-');
title(['色噪声下',num2str(Array_Num),'线阵估计',num2str(num),'个信源']);
xlabel('不同信噪比（dB）');
ylabel('正确检测概率(%)');
axis([-30 10 0 1]);
legend('AIC','MDL','GDE','BIC');

figure(2);
xx=-30:10;
plot(xx,Pd_BIC,'>-',xx,Pd_BIC_qita10,'rs-',xx,Pd_BIC_qita100,'o-');
title(['色噪声下',num2str(Array_Num),'线阵估计',num2str(num),'个信源']);
xlabel('不同信噪比（dB）');
ylabel('正确检测概率(%)');
axis([-30 10 0 1]);
legend('qita:1','qita:10','qita:100');
%%
