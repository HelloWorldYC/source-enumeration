%% 论文第五章白噪声下不同快拍数实验
clear;
clc;
tic;

f0 = 15.48e4;
fs = 62e4;
fa = 2.3e4;
fb = 2.2e4;
num=4; %信源数

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
param.errorFlag = 0;
param.errorGoal = 1e-6;
param.preserveDCAtom = 0;
% Dictionary = randn(L,param.K);
param.InitializationMethod = 'DataElements';
param.displayProgress = 0;

% L_train = 256; 
% [Dictionary_base] = construct_multidictionary(fs,L_train,fa,fb,f0,param,num_max,s_jam);
filename = strcat('./dictionaries/white/dictionary_white_sensor_', num2str(Array_Num), '.mat');
Dictionary_base = load(filename);
Dictionary_base = Dictionary_base.Dictionary_base;
%%
% 构造低通滤波器
Wp=2*pi*30;
Ws=2*pi*35;
Rp=0.5;
Rs=40;
fs1=100;
W=2*pi*fs1;
[N1,Wn]=buttord(2*Wp/W,2*Ws/W,Rp,Rs);
[b,a]=butter(N1,Wn);

Nt=200; %Monte次数
jj=0;
snr = 5;
Am=10^(snr/10);
L_circle = 30:10:300;
L_length = length(L_circle);
Pd_GDE=zeros(1,L_length);
Pd_AIC=zeros(1,L_length);
Pd_MDL=zeros(1,L_length);
Pd_IBIC=zeros(1,L_length);
Pd_ISSM=zeros(1,L_length);
Pd_MSRSE=zeros(1,L_length);
for L=L_circle
    disp(['L is ',num2str(L)]);
    jj=jj+1;
    Ns_AIC=zeros(1,Nt);
    Ns_MDL=zeros(1,Nt);
    Ns_GDE=zeros(1,Nt);
    Ns_IBIC=zeros(1,Nt);
    Ns_ISSM=zeros(1,Nt);
    Ns_MSRSE=zeros(1,Nt);
for cc=1:Nt
%     x=randn(num,L);
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
%     T1=T+sqrt(sum(T));
    [RAIC,Ns_AIC(cc)] = func_AIC(M,L,T);
    [RMDL,Ns_MDL(cc)] = func_MDL(M,L,T);
    [GDE,Ns_GDE(cc)] = func_GDE(M,L,R);
    [RBIC,Ns_IBIC(cc)] = func_IBIC(1/(M*L),M,L,R);
    [ISSM,Ns_ISSM(cc)]=func_ISSM(X);
    [MSRSE,Ns_MSRSE(cc)] = func_MSRSE(L,Dictionary_base,num_max,X,param.L);

end

Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt;
Pd_MDL(jj)=length(find(Ns_MDL==num))./Nt;
Pd_AIC(jj)=length(find(Ns_AIC==num))./Nt;
Pd_IBIC(jj)=length(find(Ns_IBIC==num))./Nt;
Pd_ISSM(jj)=length(find(Ns_ISSM==num))./Nt;
Pd_MSRSE(jj)=length(find(Ns_MSRSE==num))./Nt;

end
 %%
%  plot(L_circle,Pd_GDE,'>-',L_circle,Pd_MDL,'rs-',L_circle,Pd_AIC,'b*-',L_circle,Pd_IBIC,'r*-',L_circle,Pd_MIC,'o-',...
%      L_circle,Pd_MSTDC,'b^-',L_circle,Pd_ISSM,'gs-',L_circle,Pd_LDFCM,'rv-',L_circle,Pd_MSRSE,'ms-');
plot(L_circle,Pd_AIC,'g*-',L_circle,Pd_MDL,L_circle,Pd_IBIC,'bp-',...
    'rs-',L_circle,Pd_GDE,'m>-',L_circle,Pd_ISSM,'rd-',L_circle,Pd_MSRSE,'rs-');
xlabel('快拍数');
ylabel('正确检测概率');
axis([min(L_circle) max(L_circle) 0 1]);
legend('AIC','MDL','NBIC','GDE','ISSM','本文算法');
toc;
