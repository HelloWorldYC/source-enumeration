clear all
clc;

fs=10000;
T=1/fs;
t=0:T:1;
N=length(t);
%-----------------------------------------生成低频随机过程
fa=50;                                 %低频过程a(t)截止频率
wa=2*pi*fa/fs;                         
fb=55;                                %低频过程b(t)截止频率
wb=2*pi*fb/fs;                         
f0=2000;                                  %中心频率
at=lowfrequency(N,wa);                    %低频随机过程
figure
subplot(2,1,1),plot(t,at)
title('低频过程a'),xlabel('t'),ylabel('b(t)')
bt=lowfrequency(N,wb);
subplot(2,1,2),plot(t,bt)
title('低频过程b'),xlabel('t'),ylabel('b(t)')
%-----------------------------------------窄带随机过程及性质
X=at.*cos(2*pi*f0*t)-bt.*sin(2*pi*f0*t);  %窄带随机过程
figure,plot(t,X),title('窄带随机过程')
Rtau=xcorr(X);                            %自相关函数
tt=-N+1:N-1;
figure,plot(tt,Rtau),title('自相关函数R_x(\tau)')
Sx=fft(Rtau);                             %功率谱密度
len=length(Sx);
k=0:len-1;
w=2*pi*(k/len-1/2)*fs;
figure,plot(w/2/pi,abs(fftshift(Sx)));xlim([-2300,2300])
title('功率谱密度S_x(\omega)'),xlabel('f/Hz')
