% 生成窄带随机信号
% 输入：fs：采样率；fa：低频过程a(t)截止频率；fb：低频过程b(t)截止频率；f0：中心频率
% 输出：t：时序；at：低频过程a(t)；bt：低频过程b(t)；X：窄带随机过程；
function [t,at,bt,X]=narrow_signal(fs,L,fa,fb,f0)

T=1/fs;
t=T:T:L*T;
N=length(t);
%-----------------------------------------生成低频随机过程
wa=2*fa/fs;                         
wb=2*fb/fs;                         
at=lowfrequency(N,wa);                    %低频随机过程
bt=lowfrequency(N,wb);
%-----------------------------------------窄带随机过程及性质
X=at.*cos(2*pi*f0*t)-bt.*sin(2*pi*f0*t);  %窄带随机过程
% Rtau=xcorr(X);                            %自相关函数
% tt=-N+1:N-1;
% figure,plot(tt,Rtau),title('自相关函数R_x(\tau)')
% Sx=fft(Rtau);                             %功率谱密度
% len=length(Sx);
% k=0:len-1;
% w=2*pi*(k/len-1/2)*fs;
% figure,plot(w/2/pi,abs(fftshift(Sx)));xlim([-2300,2300])
% title('功率谱密度S_x(\omega)'),xlabel('f/Hz')

end