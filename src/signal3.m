clear all
clc;

fs=10000;
T=1/fs;
t=0:T:1;
N=length(t);
%-----------------------------------------���ɵ�Ƶ�������
fa=50;                                 %��Ƶ����a(t)��ֹƵ��
wa=2*pi*fa/fs;                         
fb=55;                                %��Ƶ����b(t)��ֹƵ��
wb=2*pi*fb/fs;                         
f0=2000;                                  %����Ƶ��
at=lowfrequency(N,wa);                    %��Ƶ�������
figure
subplot(2,1,1),plot(t,at)
title('��Ƶ����a'),xlabel('t'),ylabel('b(t)')
bt=lowfrequency(N,wb);
subplot(2,1,2),plot(t,bt)
title('��Ƶ����b'),xlabel('t'),ylabel('b(t)')
%-----------------------------------------խ��������̼�����
X=at.*cos(2*pi*f0*t)-bt.*sin(2*pi*f0*t);  %խ���������
figure,plot(t,X),title('խ���������')
Rtau=xcorr(X);                            %����غ���
tt=-N+1:N-1;
figure,plot(tt,Rtau),title('����غ���R_x(\tau)')
Sx=fft(Rtau);                             %�������ܶ�
len=length(Sx);
k=0:len-1;
w=2*pi*(k/len-1/2)*fs;
figure,plot(w/2/pi,abs(fftshift(Sx)));xlim([-2300,2300])
title('�������ܶ�S_x(\omega)'),xlabel('f/Hz')
