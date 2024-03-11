%% 信号频谱分析

clear;
clc;

j = sqrt(-1);
BW=1.5e6;
fs = 62e6;
T = 1/fs;
f0 = 0;
x1 = xlsread('E:\MYC\main\labdata\EWS_20_waveform','D101:D1100');
N = length(x1);
t = (0:N-1)*T;
X1 = fftshift(fft(x1));
f1=(0:N-1)/N*fs-fs/2;
figure;
plot(x1);
figure;
plot(f1,abs(X1));
figure;
hist(x1,N);

x2=round(50*randn(1,N));
figure;
plot(x2);
figure;
plot(abs(fft(x2)));
figure;
hist(x2,N);

x3 = normrnd(0,1,1,N);
figure;
plot(x3);
figure;
plot(abs(fft(x3)));
figure;
hist(x3,N);

% % 设计带通滤波器
fp1=[];
mu = 0;
sigma = 5;
x4 = round(50*normrnd(mu,sigma,1,N));
% x5 = round(50*normrnd(mu,sigma,1,N));
nt = x4.*cos(2*pi*f0*t)-x4.*sin(2*pi*f0*t);
fn=(0:N-1)/N*fs-fs/2;
Xn = fft(nt);
figure;
plot(nt);
figure;
plot(abs(Xn));
figure;
hist(nt,N);

mu = 0;
sigma = 5;
x6 = round(50*normrnd(mu,sigma,N,1));
X6 = fftshift(fft(x6));
figure;
hist(x6,N);
figure;
plot(abs(X6));
figure;
plot(x6);

fa = 3400;
fb = 3400;
[t,at,bt,x7]=narrow_signal(fs,N,fa,fb,f0);
X7 = fftshift(fft(x7));
figure;
plot(x7);
figure;
plot(abs(X7));
figure;
hist(x7,N);

x8 = normrnd(mu,sigma,1,N);
x9=real(x8.*(exp(j*2*pi*f0*t)));
X9 = fftshift(fft(x9));
figure;
plot(x9);
figure;
plot(abs(X9));
figure;
hist(x9,N);