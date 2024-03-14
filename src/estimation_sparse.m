%%
clear;
clc;
tic;

f0 = 15e6;
fs = 60e6;
fa = 10e6;
fb = 10e6;
L=512;
num =3; %信源�?
% M=1; 
T = 1/fs;
t = (0:L-1)*T;
mu = 0;
sigma = 1;
ii = sqrt(-1);
level = 10;
unit = eye(level);
alpha = 1200;
tau = 0;

% �?疏表示参�?
param.L = 7;
param.K = 30;
param.numIteration = 50;
param.errorFlag = 0;
param.preserveDCAtom = 0;
% Dictionary = randn(L,param.K);
param.InitializationMethod = 'DataElements';
param.displayProgress = 0;
% disp('Starting to  train the dictionary');
% [Dictionary,output] = KSVD(Y,param);
% % disp(['The KSVD algorithm retrived ',num2str(output.ratio(end)),' atoms from the original dictionary']);
% coef = full(output.CoefMatrix);
% [M,L] = size(coef);
% R = coef*coef'/L;
% [BIC,Ns_BIC(cc)] = func_BIC(1/(M*L),M,L,R);

[t1,at1,bt1,x1]=narrow_signal(fs,L,fa,fb,f0);
[t2,at2,bt2,x2]=narrow_signal(fs,L,fa,fb,f0);
[t3,at3,bt3,x3]=narrow_signal(fs,L,fa,fb,f0);
x = x1+x2+x3;
x = x';
[Yl,Yh,Yscale] = dtwavexfm(x,level,'near_sym_b','qshift_b');%Yh中从高频到低�?
for j = 1:level
    z(:,j) = dtwaveifm(Yl,Yh,'near_sym_b','qshift_b',unit(j,:));
end
% [z, u_hat, omega] = VMD(x, alpha, tau, level, false, 1, 1e-7);
train_signal = [z x];
train_signal = train_signal';
[Dictionary_base,output_base] = KSVD(train_signal,param);

% 构�?�低通滤波器
% Wp=2*pi*20;
% Ws=2*pi*40;
% Rp=0.5;
% Rs=40;
% fs=100;
% W=2*pi*fs;
% [N1,Wn]=buttord(2*Wp/W,2*Ws/W,Rp,Rs);
% [b,a]=butter(N1,Wn);

Nt=50; %Monte次数
jj=0;
snr=-20:20;
Pd_GDE=zeros(1,length(snr));
Pd_AIC=zeros(1,length(snr));
Pd_MDL=zeros(1,length(snr));
Pd_NBIC=zeros(1,length(snr));
Pd_MIC=zeros(1,length(snr));
Pd_MSTDC=zeros(1,length(snr));
Pd_GDE_sparse=zeros(1,length(snr));
Pd_BIC_sparse=zeros(1,length(snr));
for SNR=snr 
    disp(['SNR is ',num2str(SNR)]);
    source_power=10.^(SNR./10);
    source_amplitude = sqrt(source_power)*ones(1,num);    % 信源标准�?
    jj=jj+1;
    Ns_AIC=zeros(1,Nt);
    Ns_MDL=zeros(1,Nt);
    Ns_GDE=zeros(1,Nt);
    Ns_NBIC=zeros(1,Nt);
    Ns_MIC=zeros(1,Nt);
    Ns_MSTDC=zeros(1,Nt);
    Ns_GDE_sparse=zeros(1,Nt);
    Ns_BIC_sparse=zeros(1,Nt);

parfor cc=1:Nt
    z = zeros(L,level);
    [t1,at1,bt1,x1]=narrow_signal(fs,L,fa,fb,f0);
    [t2,at2,bt2,x2]=narrow_signal(fs,L,fa,fb,f0);
    [t3,at3,bt3,x3]=narrow_signal(fs,L,fa,fb,f0);
    x = [x1;x2;x3];
    source_wave = sqrt(0.5)*x;
    st=source_amplitude*source_wave;
    nt=sqrt(0.5)*randn(1,L);
%     color_noise=sqrt(0.5)*filter(b,a,nt);        %滤波产生高斯色噪�?
%     signal=st+nt;
    signal = awgn(st,SNR,'measured');
    signal = signal';

    [Yl,Yh,Yscale] = dtwavexfm(signal,level,'near_sym_b','qshift_b');%Yh中从高频到低�?
    for j = 1:level
        z(:,j) = dtwaveifm(Yl,Yh,'near_sym_b','qshift_b',unit(j,:));
    end
% 	[z, u_hat, omega] = VMD(signal, alpha, tau, level, false, 1, 1e-7);
%     z=emd(signal);
    Y = [z signal];
    Y = Y';
    [M1,~]=size(Y);
    R1 = Y*Y'/L;
    
    [~,v]=svd(R1);
    T=diag(v);
    T1=T+sqrt(sum(T));
    [AIC,Ns_AIC(cc)] = func_AIC(M1,L,T1);
    [MDL,Ns_MDL(cc)] = func_MDL(M1,L,T1);
    [GDE,Ns_GDE(cc)] = func_GDE(M1,L,R1);
    [NBIC,Ns_NBIC(cc)] = func_NBIC(1/(M1*L),M1,L,R1);
    [MIC,Ns_MIC(cc)] = func_MIC(Y,M1);
    [MSTDC,Ns_MSTDC(cc)] = func_MSTDC(Y,M1);
    
    [Dictionary,output] = KSVD(Y,param);
%     coef = full(output.CoefMatrix);
    coef = OMP(Dictionary,Y,param.L);
    resignal = Dictionary_base*coef;
    [M2,~]=size(resignal);
    R2 = resignal*resignal'/L;
    
    [GDE_sparse,Ns_GDE_sparse(cc)] = func_GDE(M2,L,R2);
    [BIC_sparse,Ns_BIC_sparse(cc)] = func_NBIC(1/(M2*L),M2,L,R2);
end

Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt;
Pd_MDL(jj)=length(find(Ns_MDL==num))./Nt;
Pd_AIC(jj)=length(find(Ns_AIC==num))./Nt;
Pd_NBIC(jj)=length(find(Ns_NBIC==num))./Nt;
Pd_MIC(jj)=length(find(Ns_MIC==num))./Nt;
Pd_MSTDC(jj)=length(find(Ns_MSTDC==num))./Nt;

Pd_GDE_sparse(jj)=length(find(Ns_GDE_sparse==num))./Nt;
Pd_BIC_sparse(jj)=length(find(Ns_BIC_sparse==num))./Nt;

end
 %%

plot(snr,Pd_GDE_sparse,'>-',snr,Pd_BIC_sparse,'rs-',snr,Pd_GDE,'o-',snr,Pd_NBIC,'b^-',snr,Pd_MIC,'gs-');
title(['单�?�道',num2str(num),'个信�?']);
xlabel('不同信噪比（dB�?');
ylabel('正确�?测概�?(%)');
axis([min(snr) max(snr) 0 1]);
legend('GDEsparse','BICsparse','GDE','BIC','MIC');
toc;
%%
% 
% Rtau = xcorr(x1);
% tt=-L+1:L-1;
% figure,plot(tt,Rtau),title('自相关函数R_x(\tau)')
% Sx=fft(Rtau);                             %功率谱密�?
% len=length(Sx);
% k=0:len-1;
% w=2*pi*(k/len-1/2)*fs;
% figure,plot(w/2/pi,abs(fftshift(Sx)));
% title('功率谱密度S_x(\omega)'),xlabel('f/Hz')