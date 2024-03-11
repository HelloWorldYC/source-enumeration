% ����խ������ź�
% ���룺fs�������ʣ�fa����Ƶ����a(t)��ֹƵ�ʣ�fb����Ƶ����b(t)��ֹƵ�ʣ�f0������Ƶ��
% �����t��ʱ��at����Ƶ����a(t)��bt����Ƶ����b(t)��X��խ��������̣�
function [t,at,bt,X]=narrow_signal(fs,L,fa,fb,f0)

T=1/fs;
t=T:T:L*T;
N=length(t);
%-----------------------------------------���ɵ�Ƶ�������
wa=2*fa/fs;                         
wb=2*fb/fs;                         
at=lowfrequency(N,wa);                    %��Ƶ�������
bt=lowfrequency(N,wb);
%-----------------------------------------խ��������̼�����
X=at.*cos(2*pi*f0*t)-bt.*sin(2*pi*f0*t);  %խ���������
% Rtau=xcorr(X);                            %����غ���
% tt=-N+1:N-1;
% figure,plot(tt,Rtau),title('����غ���R_x(\tau)')
% Sx=fft(Rtau);                             %�������ܶ�
% len=length(Sx);
% k=0:len-1;
% w=2*pi*(k/len-1/2)*fs;
% figure,plot(w/2/pi,abs(fftshift(Sx)));xlim([-2300,2300])
% title('�������ܶ�S_x(\omega)'),xlabel('f/Hz')

end