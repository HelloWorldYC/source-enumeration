function [m_t] = lowfrequency(N,omega)
%lowfrequency ������Ƶ��˹����ź�
%   t-ʱ��,w-Ƶ�ʣ�m_t -��Ƶ����
% x=wgn(1,N,5);         %������˹������
x=sqrt(5)*randn(1,N);   %��wgn������ʵ��һ����
[b,a]=butter(10,omega,'low');
m_t=filter(b,a,x);
end