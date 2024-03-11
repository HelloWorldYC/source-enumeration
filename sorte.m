clc
clear all
close all

Ns =200;
L=Ns;
num=3; %信源数
%x=randn(num,Ns);


Array_Num=7 ;% 阵元数
d=0.5; %圆阵半径
lamda=1; %波长
kk=6;
M=Array_Num;
% 入射角
theta_jam1=-30;
theta_jam2=45;
theta_jam3=60;
% theta_jam1=35;
% theta_jam2=80;
% theta_jam3=70;
% theta_jam4=70;
% theta_jam5=88;
% degrad=pi/180;
%方位角
alfa_jam1=0;
alfa_jam2=0;
alfa_jam3=90;
alfa_jam4=140;
alfa_jam5=190;

% 构造阵源方向向量
% m=[0:Array_Num-1]';
% A1=exp(-2*pi*j*(d/lamda)*m*sin(theta1*degrad));
% A2=exp(-2*pi*j*(d/lamda)*m*sin(theta2*degrad));
% A3=exp(-2*pi*j*(d/lamda)*m*sin(theta3*degrad));
% A4=exp(-2*pi*j*(d/lamda)*m*sin(theta4*degrad));
% A5=exp(-2*pi*j*(d/lamda)*m*sin(theta5*degrad));
% A=[A1,A2,A3,A4,A5];%方向矩阵；
s_jam1=array_form(Array_Num,d,lamda,theta_jam1,alfa_jam1,kk);
s_jam2=array_form(Array_Num,d,lamda,theta_jam2,alfa_jam2,kk);
s_jam3=array_form(Array_Num,d,lamda,theta_jam3,alfa_jam3,kk);
% s_jam4=array_form(Array_Num,d,lamda,theta_jam4,alfa_jam4,kk);
% s_jam5=array_form(Array_Num,d,lamda,theta_jam5,alfa_jam5,kk);
A=[s_jam1;s_jam2;s_jam3];%方向矩阵；
A = A';

Nt=1000;
jj=0;

for SNR= -20:20

    jj=jj+1;
    Am=10^(SNR/10);

for cc=1:Nt

  x=randn(num,Ns);
  signal=Am.*x; 
  A1=A*signal;%模拟天线接收到的信号
 noise=randn(M,L)+j*randn(M,L);


%  X = awgn(A1,SNR,1);%白噪声模型
X=A1+noise;
R=X*X'/L; %信号协方差
[u,v]=svd(R);

D=diag(v);
D=sort(D,'descend');
lambda_1=zeros(1,M-1);
for i=1:M-1
    lambda_1(i)=D(i)-D(i+1);
end
var_k=zeros(1,M-1);
for i=1:M-1
    var_k(i)=var(lambda_1(i:M-1));
end
for i=1:M-2
    if(var_k(i)~=0)
        SORTE(i)=var_k(i+1)/var_k(i);
    else
        SORTE(i)=1e10;
    end
end
 [aa num_1]=sort(SORTE(1:M-3));
 SORTE(cc)=num_1(1);%信源数
end
Pd_SORTE(jj)=length(find(SORTE==num))./Nt;
end
xx=-20:20;
plot(xx,Pd_SORTE);

