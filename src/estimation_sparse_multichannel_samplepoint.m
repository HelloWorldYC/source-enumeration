%% 利用稀疏矩阵生成字典，再估计
clear;
clc;
tic;

f0 = 15.48e4;
fs = 62e4;
fa = 2.3e4;
fb = 2.2e4;
num=2; %信源数

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
Dictionary_base = load('dictionary_color.mat');
Dictionary_base = Dictionary_base.Dictionary_base;
%%
Wp=2*pi*30;
Ws=2*pi*40;
Rp=0.5;
Rs=40;
fs1=100;
W=2*pi*fs1;
[N1,Wn]=buttord(2*Wp/W,2*Ws/W,Rp,Rs);
[b,a]=butter(N1,Wn);

Nt=100; %Monte次数
jj=0;
snr = 10;
Am=10^(snr/10);
L_circle = 30:10:400;
L_length = length(L_circle);
Pd_GDE=zeros(1,L_length);
Pd_AIC=zeros(1,L_length);
Pd_MDL=zeros(1,L_length);
Pd_IBIC=zeros(1,L_length);
Pd_MIC=zeros(1,L_length);
Pd_MSTDC=zeros(1,L_length);
Pd_ISSM=zeros(1,L_length);
Pd_LDFCM=zeros(1,L_length);
Pd_MSRSE=zeros(1,L_length);
% Pd_GDE_sparse=zeros(1,L_length);
% Pd_BIC_sparse=zeros(1,L_length);
% Pd_MIC_sparse=zeros(1,L_length);
for L=L_circle
    disp(['L is ',num2str(L)]);
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
%     Ns_GDE_sparse=zeros(1,Nt);
%     Ns_BIC_sparse=zeros(1,Nt);
%     Ns_MIC_sparse=zeros(1,Nt);
parfor cc=1:Nt
%     x=randn(num,L);
    x1 = zeros(num,L);
    for i=1:num
        [t1,at1,bt1,x1(i,:)]=narrow_signal(fs,L,fa,fb,f0);
    end
%     [code1,x1] = BPSK_generate(500,1000,2000,128,1);
%     [code2,x2] = BPSK_generate(500,1000,2000,128,1);
%     [code3,x3] = BPSK_generate(500,1000,2000,128,1);
    signal=Am*x1;
    A1=A*signal; 
    X=awgn(A1,snr,'measured');
%     noise=randn(M,L);
%     color_noise=filter(b,a,noise);        %滤波产生高斯色噪声
%     X=A1+color_noise;

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
%     [Dictionary,output] = KSVD(X,param);
%     coef = OMP(Dictionary_base,X,param.L);
%     resignal = Dictionary_base*coef;
%     [M2,~]=size(resignal);
%     R2 = resignal*resignal'/L;
%     
%     [GDE_sparse,Ns_GDE_sparse(cc)] = func_GDE(M2,L,R2);
%     [BIC_sparse,Ns_BIC_sparse(cc)] = func_BIC(1/(M2*L),M2,L,R2);
%     [MIC_sparse,Ns_MIC_sparse(cc)] = func_MIC(resignal,M2,L);

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
% Pd_GDE_sparse(jj)=length(find(Ns_GDE_sparse==num))./Nt;
% Pd_BIC_sparse(jj)=length(find(Ns_BIC_sparse==num))./Nt;
% Pd_MIC_sparse(jj)=length(find(Ns_MIC_sparse==num))./Nt;

end
 %%
 plot(L_circle,Pd_GDE,'>-',L_circle,Pd_MDL,'rs-',L_circle,Pd_AIC,'b*-',L_circle,Pd_IBIC,'r*-',L_circle,Pd_MIC,'o-',...
     L_circle,Pd_MSTDC,'b^-',L_circle,Pd_ISSM,'gs-',L_circle,Pd_LDFCM,'rv-',L_circle,Pd_MSRSE,'ms-');
% plot(L_circle,Pd_AIC,'g*-',L_circle,Pd_MDL,'bp-',L_circle,Pd_GDE,'m>-',...
%      L_circle,Pd_MSTDC,'go-',L_circle,Pd_IBIC,'b^-',L_circle,Pd_ISSM,'md-',L_circle,Pd_MSRSE,'rs-');
% plot(L_circle,Pd_GDE_sparse,'>-',L_circle,Pd_BIC_sparse,'rs-',L_circle,Pd_MIC_sparse,'r*-',L_circle,Pd_GDE,'o-',L_circle,Pd_BIC,'b^-',L_circle,Pd_MIC,'gs-');
% plot(xx,Pd_AIC,'>-',xx,Pd_MDL,'rs-',xx,Pd_GDE,'o-',xx,Pd_BIC,'b^-',xx,Pd_MIC,'c*-',xx,Pd_MSTDC,'r*-');
% title(['白/色噪声下',num2str(Array_Num),'线阵估计',num2str(num),'个信源']);
xlabel('Number of Samples');
ylabel('Probability of Correct Detection');
axis([min(L_circle) max(L_circle) 0 1]);
legend('GDE','MDL','AIC','IBIC','MIC','MSTDC','ISSM','LDFCM','proposed');
% legend('AIC','MDL','GDE','MSTDC','NBIC','ISSM','proposed');
toc;
