clear ;
close all;
clc;%�������

T=100e-6;%����ʱ��
fs=300e6;%����Ƶ��
N=T*fs;%��������
detlf=20e6;%�˲�����ֹƵ��
f1=100e6;%�����ź�����Ƶ��
m=0.5;%���ƶ�
kfm=5e6;%��Ƶб��
kpm=5;%����б��
M=100;%���۴���
p=fft(fir1(N-1,detlf/fs*2));%�˲���Ƶ��
s=0;

for i=1:100
    xn=ifft(fft(random('Normal',0,1,1,N)).*p);%��˹������ͨ���˲���
    j=abs(fft(xn));
    s=s+j;
end

s=s/M;
j=s;

figure(1)
t=0:1/fs:T-1/fs;
plot(t*1e6,xn);
xlabel('us');
title('��Ƶ����ʱ����');

figure(2)
f=(0:N-1)*fs/N;
plot(f*1e-6,20*log10(j.^2/max(j.^2)));%��Ƶ����������
axis([-1 22 -8 0]);
xlabel('MHZ');
title('��Ƶ����������');

n=1:N;
zn=(1+m*cos(2*pi*xn)).*cos(2*pi*f1/fs*n);%�����������ű��ʽ

figure(3)
plot(t*1e6,zn);
title('������������ʱ����');
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
plot(f*1e-6,20*log10(j.^2/max(j.^2)));%�����������Ź�����
title('�����������Ź�����');
xlabel('MHZ');
axis([90 110 -200 0]);

sum(1)=0;
for i=1:N-1;
    sum(i+1)=xn(i)+sum(i);
end

xn=sum/fs;
wn=cos((2*pi*f1*t+2*pi*kfm*xn));%������Ƶ���ű��ʽ

figure(5)
plot(t*1e6,wn);
title('������Ƶ����ʱ����');
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
plot(f*1e-6,20*log10(j.^2/max(j.^2)));%������Ƶ���Ź�����
axis([50 150 -150 0])
xlabel('MHZ');
title('������Ƶ���Ź�����');

sum(1)=0;
for i=1:N-1;
    sum(i+1)=xn(i)+sum(i);
end

xn=sum/fs;
on=cos(2*pi*f1*t+kpm*xn);%����������ű��ʽ

figure(7)
plot(t*1e6,on);
title('�����������ʱ����');
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
plot(f*1e-6,20*log10(j.^2/max(j.^2)));%����������Ź�����
axis([50 150 -70 0])
xlabel('MHZ');
title('����������Ź�����');