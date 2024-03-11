%% PSO优化VMD测试

clear;
clc;
close all;
tic;

jj = sqrt(-1);
fs = 1000;
ts = 1/fs;
L = 500;
t = (0:L-1)*ts;
num = 3;
d=0.5; %线阵半径
lamda=1; %波长
kk=6;    %线阵
% 入射角
theta_jam1=10;
theta_jam2=40;
theta_jam3=70;
%方位角
alfa_jam1=10;
alfa_jam2=50;
alfa_jam3=90;

s_jam1=array_form(8,d,lamda,theta_jam1,alfa_jam1,kk);
s_jam2=array_form(8,d,lamda,theta_jam2,alfa_jam2,kk);
s_jam3=array_form(8,d,lamda,theta_jam3,alfa_jam3,kk);
A=[s_jam1;s_jam2;s_jam3];%阵列流型；
A=A';

snr = 20;
Am=10^(snr/10);
x = randn(num,L);
% x=normrnd(3,1.5,[num,L]);
signal=Am*x;
A1=A*signal;
A2 = A1(1,:);
X = awgn(A2,snr,'measured');

%
% snr = 10;
% SNR = ones(num,1)*snr;
% source_power = (10.^(SNR./10));
% source_amplitude = sqrt(source_power)*ones(1,L); %信源标准差
% source_wave = sqrt(0.5)*(randn(L,num)+jj*randn(L,num));
% st = source_amplitude.*source_wave';
% nt = sqrt(0.5)*(randn(1,L)+jj*randn(1,L));
% X = A*st+nt;

% % 线性调频信号
% f0=100;
% f1=80;
% f2=120;
% c=1500;
% lambda=c/f0;
% d=lambda/2;
% SNR=15;
% b=pi/180;
% theat1=30*b;    %入射信号波束角1
% theat2=0*b;     %入射信号波束角2
% n=ts:ts:L*ts;
% theat=[theat1 theat2];
% s1=chirp(n,80,1,120);%生成线性调频信号1
% sa=fft(s1,2048);
% s2=chirp(n+0.100,40,1,200);%生成线性调频信号2


%%
figure(1);
subplot(4,1,1);
plot(x(1,:));
subplot(4,1,2);
plot(x(2,:));
subplot(4,1,3);
plot(x(3,:));
subplot(4,1,4);
plot(X);

%%
% num_range = [4,10];
% alpha_range = [0,3000];
% range = [num_range;alpha_range];
% max_v = 0.2*(range(:,2)-range(:,1));
% tau = 0;
% pso_Trelea_vectorized('PSO_VMD_fun',2,max_v,range);
% [Y, u_hat, omega] = VMD(X, alpha, tau, num_mode, false, 1, 1e-7);
alpha = 1200;
tau = 0;
num_mode = 6;
[Y, u_hat, omega] = VMD(X, alpha, tau, num_mode, false, 2, 1e-7);

%%
figure(3);
for i =1:num_mode
   subplot(num_mode,1,i);
   plot(Y(i,:));
end
corrcoef(Y');
figure(4);
for i =1:num_mode
%    temp=abs(fft(Y(i,:)));
   subplot(num_mode,1,i);
%    plot((1:L/2),temp(1:L/2));
    plot(abs(fft(Y(i,:))));
end

%%
[M,L]=size(Y);
% Y1=zeros(num_mode/2,L);
% for j=1:num_mode/2
%     Y1(j,:) = Y(j,:)+Y(num_mode-j+1,:);
% end
% Y2=Y1;
% [M,L] = size(Y2);
Y2=Y;

figure(5);
for i =1:M
   subplot(M,1,i);
   plot(Y2(i,:));
end
%%
% 自相关系数矩阵
% y_mean = mean(Y2,2);
% for v=1:M
%    Y2(v,:)=Y2(v,:)-y_mean(v); 
% end
% V=Y2*Y2'/L;
% V_diag=diag(diag(V));
% R = (V_diag^(-1/2))*V*(V_diag^(-1/2));

R = Y2*Y2'/L;
[BIC,Ns_BIC] = func_BIC(1/(M*L),M,L,R);
[GDE,Ns_GDE] = func_GDE(M,L,R);
