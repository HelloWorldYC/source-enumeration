%% 尝试用稀疏表示增强现有信源数估计算法
clear;
clc;

f0 = 15.48e4;
fs = 62e4;
fa = 2.3e4;
fb = 2.2e4;
Ns = 256;
L=Ns;
t=1:L;
num=3; %信源数

Array_Num=8;% 阵元数
d=0.5; %线阵半径
lamda=1; %波长
kk=6;    %线阵
M=Array_Num;
num_max = 6;
% 入射角
theta_jam=10:15:num_max*20;
degrad=pi/180;
%方位角
alfa_jam=10:20:num_max*20;

s_jam = zeros(num_max,M);
for i=1:num_max
s_jam(i,:)=array_form(Array_Num,d,lamda,theta_jam(i),alfa_jam(i),kk);
end
A=s_jam(1:num,:);%方向矩阵；
A=A';

% 稀疏表示参数
param.L = 3;
param.K = 45;
param.numIteration = 50;
param.errorFlag = 1;
param.errorGoal = 1e-6;
param.preserveDCAtom = 0;
% Dictionary = randn(L,param.K);
param.InitializationMethod = 'DataElements';
param.displayProgress = 0;

% [Dictionary_base] = construct_multidictionary(fs,L,fa,fb,f0,param,num_max,s_jam);
filename = strcat('./dictionaries/white/dictionary_white_sensor_', num2str(Array_Num), '.mat');
% filename = strcat('dictionary_white.mat');
Dictionary_base = load(filename);
Dictionary_base = Dictionary_base.Dictionary_base;

%%
Nt=10; %Monte次数
jj=0;
snr_min = -20;
snr_max = 20;
snr_length = snr_max-snr_min+1;
Pd_AIC=zeros(1,snr_length);
Pd_MDL=zeros(1,snr_length);
Pd_MDL_sparse=zeros(1,snr_length);
Pd_NBIC=zeros(1,snr_length);
Pd_GDE=zeros(1,snr_length);
Pd_ISSM=zeros(1,snr_length);
Pd_ISSM_sparse=zeros(1,snr_length);
Pd_NBIC_sparse=zeros(1,snr_length);
for SNR=snr_min:snr_max 
    disp(['SNR is ',num2str(SNR)]);
    Am=10^(SNR/10);
    jj=jj+1;
    Ns_AIC=zeros(1,Nt);
    Ns_MDL=zeros(1,Nt);
    Ns_MDL_sparse=zeros(1,Nt);
    Ns_NBIC=zeros(1,Nt);
    Ns_GDE=zeros(1,Nt);
    Ns_ISSM=zeros(1,Nt);
    Ns_NBIC_sparse=zeros(1, Nt);
    Ns_ISSM_sparse=zeros(1, Nt);
for cc=1:Nt
    x1 = zeros(num,L);
    for i=1:num
        [t1,at1,bt1,x1(i,:)]=narrow_signal(fs,L,fa,fb,f0);
    end
    signal=Am*x1;
    A1=A*signal; 
    X=awgn(A1,SNR,'measured');

    R=X*X'/L; %信号协方差
    [u,v]=svd(R);
    T=diag(v);
    [AIC,Ns_AIC(cc)] = func_AIC(M,L,T);
    [MDL,Ns_MDL(cc)] = func_MDL(M,L,T);
    [NBIC,Ns_NBIC(cc)] = func_NBIC(1/(M*L),M,L,R);
    [GDE,Ns_GDE(cc)] = func_GDE(M,L,R);
    [ISSM,Ns_ISSM(cc)]=func_ISSM(X);
    coef = OMP(Dictionary_base{num},X,param.L);
    coef = full(coef);
    resignal = Dictionary_base{num}*coef;
    reR = resignal*resignal'/L;
    [uRe,vRe]=svd(reR);
    TRe=diag(vRe);
    reR=resignal*resignal'/L; 
    [~,Ns_MDL_sparse(cc)] = func_MDL(M,L,TRe);
    [~, Ns_NBIC_sparse(cc)]=func_NBIC(1/(M*L),M,L,reR);
    [~,Ns_ISSM_sparse(cc)]=func_ISSM(X);
end

Pd_MDL(jj)=length(find(Ns_MDL==num))./Nt;
Pd_MDL_sparse(jj)=length(find(Ns_MDL_sparse==num))./Nt;
Pd_NBIC(jj)=length(find(Ns_NBIC==num))./Nt;
Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt; 
Pd_ISSM(jj)=length(find(Ns_ISSM==num))./Nt;
Pd_NBIC_sparse(jj)=length(find(Ns_NBIC_sparse==num))./Nt;
Pd_ISSM_sparse(jj)=length(find(Ns_ISSM_sparse==num))./Nt;

end

%%
rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

xx=snr_min:snr_max;

hold on;
plot(xx,Pd_MDL,'Color',rgbTriplet(1,:),'Marker','p');
plot(xx,Pd_MDL_sparse,'Color',rgbTriplet(2,:),'Marker','p');
plot(xx,Pd_NBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(xx,Pd_NBIC_sparse,'Color',rgbTriplet(4,:),'Marker','^');
plot(xx,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');
plot(xx,Pd_ISSM_sparse,'Color',rgbTriplet(6,:),'Marker','s');

box on;
grid on;
xlabel('信噪比(dB)');
ylabel('正确检测概率');
axis([snr_min snr_max 0 1]);
legend('MDL','MDL-Enhance','NBIC','NBIC-Enhance','ISSM','ISSM-Enhance','Location','southeast');
