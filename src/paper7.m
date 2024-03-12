%% 论文第三章色噪声下不同信源数实验
clear;
clc;
tic;

f0 = 15.48e4;
fs = 62e4;
fa = 2.3e4;
fb = 2.3e4;
L = 256;

Array_Num=8;% 阵元数
d=0.5; %线阵半径
lamda=1; %波长
kk=6;    %线阵
M=Array_Num;
num_max = 6;
% 入射角
theta_jam=10:15:num_max*20;
degrad=pi/180;
%方位角
alfa_jam=10:20:num_max*20;

s_jam = zeros(num_max,M);
for i=1:num_max
s_jam(i,:)=array_form(Array_Num,d,lamda,theta_jam(i),alfa_jam(i),kk);
end

%%
% 构造低通滤波器
Wp=2*pi*35;
Ws=2*pi*40;
Rp=0.5;
Rs=40;
fs1=100;
W=2*pi*fs1;
[N1,Wn]=buttord(2*Wp/W,2*Ws/W,Rp,Rs);
[b,a]=butter(N1,Wn);

Nt=5000; %Monte次数
jj=0;
snr = 10;
Am=10^(snr/10);
num_circle = 1:1:num_max;
num_length = length(num_circle);
Pd_GDE=zeros(1,num_length);
Pd_RAIC=zeros(1,num_length);
Pd_RMDL=zeros(1,num_length);
Pd_RIBIC=zeros(1,num_length);
Pd_ISSM=zeros(1,num_length);

for num=num_circle
    A=s_jam(1:num,:);%方向矩阵；
    A=A';
    disp(['num is ',num2str(num)]);
    jj=jj+1;
    Ns_RAIC=zeros(1,Nt);
    Ns_RMDL=zeros(1,Nt);
    Ns_GDE=zeros(1,Nt);
    Ns_RIBIC=zeros(1,Nt);
    Ns_ISSM=zeros(1,Nt);

for cc=1:Nt
    x1 = zeros(num,L);
    for i=1:num
        [t1,at1,bt1,x1(i,:)]=narrow_signal(fs,L,fa,fb,f0);
    end
    signal=Am*x1;
    A1=A*signal; 
    noise=randn(M,L);
    color_noise=filter(b,a,noise);        %滤波产生高斯色噪声
    X=A1+color_noise;

    R=X*X'/L; %信号协方差

    [u,v]=svd(R);
    T=diag(v);
    T1=T+sqrt(sum(T));
    [RAIC,Ns_RAIC(cc)] = func_AIC(M,L,T1);
    [RMDL,Ns_RMDL(cc)] = func_MDL(M,L,T1);
    [GDE,Ns_GDE(cc)] = func_GDE(M,L,R);
    [RBIC,Ns_RIBIC(cc)] = func_RIBIC(1/(M*L),M,L,R);
    [ISSM,Ns_ISSM(cc)]=func_ISSM(X);

end

Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt;
Pd_RMDL(jj)=length(find(Ns_RMDL==num))./Nt;
Pd_RAIC(jj)=length(find(Ns_RAIC==num))./Nt;
Pd_RIBIC(jj)=length(find(Ns_RIBIC==num))./Nt;
Pd_ISSM(jj)=length(find(Ns_ISSM==num))./Nt;

end
%%
rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

hold on;
plot(num_circle,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(num_circle,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(num_circle,Pd_RIBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(num_circle,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(num_circle,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');

box on;
xlabel('信号源数');
ylabel('正确检测概率');
axis([min(num_circle) max(num_circle) 0 1]);
legend('RAIC','RMDL','RNBIC','GDE','ISSM','Location','southwest');
toc;
