%% 利用稀疏矩阵生成字典，再估计
clear;
clc;
% tic;

f0 = 15.48e4;
fs = 62e4;
fa = 2.3e4;
fb = 2.2e4;
Ns = 256;
L=Ns;
t=1:L;
num=4; %信源数

Array_Num=26;% 阵元数
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
% filename = strcat('./dictionaries/white/dictionary_white_sensor_', num2str(Array_Num), '.mat');
% Dictionary_base = load(filename);
% Dictionary_base = Dictionary_base.Dictionary_base;
%%
Nt=200; %Monte次数
jj=0;
snr_min = -10;
snr_max = 10;
snr_length = snr_max-snr_min+1;
Pd_GDE=zeros(1,snr_length);
Pd_AIC=zeros(1,snr_length);
Pd_MDL=zeros(1,snr_length);
Pd_IBIC=zeros(1,snr_length);
Pd_MIC=zeros(1,snr_length);
Pd_MSTDC=zeros(1,snr_length);
Pd_ISSM=zeros(1,snr_length);
Pd_LDFCM=zeros(1,snr_length);
Pd_MSRSE=zeros(1,snr_length);
Pd_test=zeros(1,snr_length);
coef = cell(1,num_max);
for SNR=snr_min:snr_max 
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
    Ns_MSRSE=zeros(1,Nt);
    Ns_test=zeros(1,Nt);
parfor cc=1:Nt
%     x=randn(num,L);
    x1 = zeros(num,L);
    for i=1:num
        [t1,at1,bt1,x1(i,:)]=narrow_signal(fs,L,fa,fb,f0);
    end
    signal=Am*x1;
    A1=A*signal; 
    X=awgn(A1,SNR,'measured');
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
    [BIC,Ns_IBIC(cc)] = func_IBIC(1/(M*L),M,L,R);
    [MIC,Ns_MIC(cc)] = func_MIC(X,M,L);
    [MSTDC,Ns_MSTDC(cc)] = func_MSTDC(X,M,L);
    [ISSM,Ns_ISSM(cc)]=func_ISSM(X);
    [Ns_LDFCM(cc)]=func_LDFCM(R);
    [MSRSE,Ns_MSRSE(cc)] = func_MSRSE(L,Dictionary_base,num_max,X,param.L);
end

Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt; 
Pd_MDL(jj)=length(find(Ns_MDL==num))./Nt;
Pd_AIC(jj)=length(find(Ns_AIC==num))./Nt;
Pd_IBIC(jj)=length(find(Ns_IBIC==num))./Nt;
Pd_MIC(jj)=length(find(Ns_MIC==num))./Nt;
Pd_MSTDC(jj)=length(find(Ns_MSTDC==num))./Nt;
Pd_ISSM(jj)=length(find(Ns_ISSM==num))./Nt;
Pd_LDFCM(jj)=length(find(Ns_LDFCM==num))./Nt;
Pd_MSRSE(jj)=length(find(Ns_MSRSE==num))./Nt;

end
%%
xx=snr_min:snr_max;
plot(xx,Pd_AIC,'g*-',xx,Pd_MDL,'bp-',xx,Pd_GDE,'m>-',...
     xx,Pd_MSTDC,'go-',xx,Pd_IBIC,'b^-',xx,Pd_ISSM,'md-',xx,Pd_MSRSE,'rs-');
xlabel('SNR(dB)');
ylabel('Probability of Correct Detection');
axis([snr_min snr_max 0 1]);
legend('AIC','MDL','GDE','MSTDC','NBIC','ISSM','proposed');
% toc;