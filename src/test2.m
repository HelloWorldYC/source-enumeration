%%
clear;
clc;
tic;

f0 = 15e6;
fs = 60e6;
fa = 8e6;
fb = 8e6;
L=500;
num =3; %信源数
% M=1; 
T = 1/fs;
t = (0:L-1)*T;
mu = 0;
sigma = 1;
ii = sqrt(-1);
alpha = 1000;
tau = 0;

% 构造低通滤波器
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
Pd_BIC=zeros(1,length(snr));
Pd_MIC=zeros(1,length(snr));
Pd_MSTDC=zeros(1,length(snr));
level = 10;
unit = eye(level);
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
parfor cc=1:Nt
    [t1,at1,bt1,x1]=narrow_signal(fs,L,fa,fb,f0);
    [t2,at2,bt2,x2]=narrow_signal(fs,L,fa,fb,f0);
    [t3,at3,bt3,x3]=narrow_signal(fs,L,fa,fb,f0);
    x = [x1;x2;x3];
    source_wave = sqrt(0.5)*x;
    st=source_amplitude*source_wave;
    nt=sqrt(0.5)*randn(1,L);
%     color_noise=sqrt(0.5)*filter(b,a,nt);        %滤波产生高斯色噪声
%     signal=st+nt;
    signal = awgn(st,SNR,'measured');
    signal = signal';

%     [Yl,Yh,Yscale] = dtwavexfm(signal,level,'near_sym_b','qshift_b');%Yh中从高频到低频
%     for j = 1:level
%         z(:,j) = dtwaveifm(Yl,Yh,'near_sym_b','qshift_b',unit(j,:));
%     end
%     z=emd(signal);
    [z, u_hat, omega] = VMD(signal, alpha, tau, level, false, 1, 1e-7); % VMD从低频到高频
    Y = [signal z'];
    Y = Y';
    [M,~]=size(Y);
    R=Y*Y'/L; %信号协方差
    
    [~,v]=svd(R);
    T=diag(v);
    T1=T+sqrt(sum(T));
    [AIC,Ns_AIC(cc)] = func_AIC(M,L,T1);
    [MDL,Ns_MDL(cc)] = func_MDL(M,L,T1);
    [GDE,Ns_GDE(cc)] = func_GDE(M,L,R);
    [BIC,Ns_BIC(cc)] = func_BIC(1/(M*L),M,L,R);
    [MIC,Ns_MIC(cc)] = func_MIC(Y,M);
    [MSTDC,Ns_MSTDC(cc)] = func_MSTDC(Y,M);

end

Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt;
Pd_MDL(jj)=length(find(Ns_MDL==num))./Nt;
Pd_AIC(jj)=length(find(Ns_AIC==num))./Nt;
Pd_BIC(jj)=length(find(Ns_BIC==num))./Nt;
Pd_MIC(jj)=length(find(Ns_MIC==num))./Nt;
Pd_MSTDC(jj)=length(find(Ns_MSTDC==num))./Nt;
end
 %%

plot(snr,Pd_AIC,'>-',snr,Pd_MDL,'rs-',snr,Pd_GDE,'o-',snr,Pd_BIC,'b^-');
title(['单通道',num2str(num),'个信源']);
xlabel('不同信噪比（dB）');
ylabel('正确检测概率(%)');
axis([min(snr) max(snr) 0 1]);
legend('AIC','MDL','GDE','BIC');
toc;

%%

% Rtau = xcorr(x1);
% tt=-L+1:L-1;
% figure,plot(tt,Rtau),title('自相关函数R_x(\tau)')
% Sx=fft(Rtau);                             %功率谱密度
% len=length(Sx);
% k=0:len-1;
% w=2*pi*(k/len-1/2)*fs;
% figure,plot(w/2/pi,abs(fftshift(Sx)));
% title('功率谱密度S_x(\omega)'),xlabel('f/Hz')