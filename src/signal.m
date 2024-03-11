clear ;
close all;
clc;%清除变量

T=100e-6;%采样时间
fs=300e6;%采样频率
N=T*fs;%采样点数
detlf=20e6;%滤波器截止频率
f1=100e6;%调制信号中心频率
m=0.5;%调制度
kfm=5e6;%调频斜率
kpm=5;%调相斜率
M=100;%积累次数
p=fft(fir1(N-1,detlf/fs*2));%滤波器频谱
s=0;

for i=1:100
    xn=ifft(fft(random('Normal',0,1,1,N)).*p);%高斯白噪声通过滤波器
    j=abs(fft(xn));
    s=s+j;
end

s=s/M;
j=s;

figure(1)
t=0:1/fs:T-1/fs;
plot(t*1e6,xn);
xlabel('us');
title('射频噪声时域波形');

figure(2)
f=(0:N-1)*fs/N;
plot(f*1e-6,20*log10(j.^2/max(j.^2)));%视频噪声功率谱
axis([-1 22 -8 0]);
xlabel('MHZ');
title('射频噪声功率谱');

n=1:N;
zn=(1+m*cos(2*pi*xn)).*cos(2*pi*f1/fs*n);%噪声调幅干扰表达式

figure(3)
plot(t*1e6,zn);
title('噪声调幅干扰时域波形');
xlabel('us');

s=0;
    for i=1:100
    zn=(1+m*cos(2*pi*xn)).*cos(2*pi*f1/fs*n);
    j=abs(fft(zn));
    s=s+j;
end

s=s/M;
j=s;

figure(4)
plot(f*1e-6,20*log10(j.^2/max(j.^2)));%噪声调幅干扰功率谱
title('噪声调幅干扰功率谱');
xlabel('MHZ');
axis([90 110 -200 0]);

sum(1)=0;
for i=1:N-1;
    sum(i+1)=xn(i)+sum(i);
end

xn=sum/fs;
wn=cos((2*pi*f1*t+2*pi*kfm*xn));%噪声调频干扰表达式

figure(5)
plot(t*1e6,wn);
title('噪声调频干扰时域波形');
xlabel('us');

s=0;
for i=1:10
    xn=ifft(fft(random('Normal',0,1,1,N)).*p);
    sum(1)=0;
    for i=1:N-1;
        sum(i+1)=xn(i)+sum(i);
    end
    xn=sum/fs;
    wn=cos((2*pi*f1*t+2*pi*kfm*xn));
    j=abs(fft(wn));
    s=s+j;
end

s=s/M;
j=s;

figure(6)
plot(f*1e-6,20*log10(j.^2/max(j.^2)));%噪声调频干扰功率谱
axis([50 150 -150 0])
xlabel('MHZ');
title('噪声调频干扰功率谱');

sum(1)=0;
for i=1:N-1;
    sum(i+1)=xn(i)+sum(i);
end

xn=sum/fs;
on=cos(2*pi*f1*t+kpm*xn);%噪声调相干扰表达式

figure(7)
plot(t*1e6,on);
title('噪声调相干扰时域波形');
xlabel('us');

s=0;
for i=1:100
    xn = ifft(fft(random('Normal',0,1,1,N)).*p);
    on=cos(2*pi*f1*t+kpm*xn);
    j=abs(fft(on));
    s=s+j;
end

s=s/M;
j=s;

figure(8)
plot(f*1e-6,20*log10(j.^2/max(j.^2)));%噪声调相干扰功率谱
axis([50 150 -70 0])
xlabel('MHZ');
title('噪声调相干扰功率谱');