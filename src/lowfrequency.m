function [m_t] = lowfrequency(N,omega)
%lowfrequency 产生低频高斯随机信号
%   t-时间,w-频率，m_t -低频噪声
% x=wgn(1,N,5);         %产生高斯白噪声
x=sqrt(5)*randn(1,N);   %跟wgn生成其实是一样的
[b,a]=butter(10,omega,'low');
m_t=filter(b,a,x);
end