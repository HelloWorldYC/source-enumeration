%% 论文第五章白噪声下不同阵元数实验
clear;
clc;
% tic;

f0 = 15.48e4;
fs = 62e4;
fa = 2.3e4;
fb = 2.2e4;
Ns = 256;
L=Ns;
t=1:L;
num=4; %信源数

% 构造低通滤波器
Wp=2*pi*30;
Ws=2*pi*35;
Rp=0.5;
Rs=40;
fs1=100;
W=2*pi*fs1;
[N1,Wn]=buttord(2*Wp/W,2*Ws/W,Rp,Rs);
[b,a]=butter(N1,Wn);

sensor_min = 8;
sensor_max = 30;
d=0.5; %线阵半径
lamda=1; %波长
kk=6;    %线阵
num_max = 6;
% 入射角
theta_jam=10:15:num_max*20;
degrad=pi/180;
%方位角
alfa_jam=10:20:num_max*20;

% 稀疏表示参数
param.L = 3;
param.K = 45;
param.numIteration = 50;
param.errorFlag = 1; 
param.errorGoal = 1e-8;
param.preserveDCAtom = 0;
% Dictionary = randn(L,param.K);
param.InitializationMethod = 'DataElements';
param.displayProgress = 0;

for Array_Num=sensor_min:sensor_max % 阵元数
disp(['ArrayNum is ',num2str(Array_Num)]);

M=Array_Num;

s_jam = zeros(num_max,M);
for i=1:num_max
s_jam(i,:)=array_form(Array_Num,d,lamda,theta_jam(i),alfa_jam(i),kk);
end
A=s_jam(1:num,:);%方向矩阵；
A=A';

end

%%
Nt=5000; %Monte次数
jj=0;
SNR = 3;
sensor_length = sensor_max-sensor_min + 1;
Pd_GDE=zeros(1,sensor_length);
Pd_AIC=zeros(1,sensor_length);
Pd_MDL=zeros(1,sensor_length);
Pd_IBIC=zeros(1,sensor_length);
Pd_ISSM=zeros(1,sensor_length);
Pd_MSRSE=zeros(1,sensor_length);

coef = cell(1,num_max);
for sensor=sensor_min:1:sensor_max
    filename = strcat('./dictionaries/white/dictionary_white_sensor_', num2str(sensor), '.mat');
    Dictionary_base = load(filename);
    Dictionary_base = Dictionary_base.Dictionary_base;
    disp(['ArrayNum is ',num2str(sensor)]);
    
    s_jam = zeros(num_max,sensor);
    for i=1:num_max
        s_jam(i,:)=array_form(sensor,d,lamda,theta_jam(i),alfa_jam(i),kk);
    end
    A=s_jam(1:num,:);%方向矩阵；
    
    A=A';
    
    Am=10^(SNR/10);
    jj=jj+1;
    Ns_AIC=zeros(1,Nt);
    Ns_MDL=zeros(1,Nt);
    Ns_GDE=zeros(1,Nt);
    Ns_IBIC=zeros(1,Nt);
    Ns_ISSM=zeros(1,Nt);
    Ns_MSRSE=zeros(1,Nt);
for cc=1:Nt
    x1 = zeros(num,L);
    for i=1:num
        [t1,at1,bt1,x1(i,:)]=narrow_signal(fs,L,fa,fb,f0);
    end
    signal=Am*x1;
    A1=A*signal; 
    X=awgn(A1,SNR,'measured');

    R=X*X'/L; %信号协方差

    [u,v]=svd(R);
    T=diag(v);
    [RAIC,Ns_AIC(cc)] = func_AIC(sensor,L,T);
    [RMDL,Ns_MDL(cc)] = func_MDL(sensor,L,T);
    [GDE,Ns_GDE(cc)] = func_GDE(sensor,L,R);
    [RBIC,Ns_IBIC(cc)] = func_IBIC(1/(sensor*L),sensor,L,R);
    [ISSM,Ns_ISSM(cc)]=func_ISSM(X);
    [MSRSE,Ns_MSRSE(cc)] = func_MSRSE(L,Dictionary_base,num_max,X,param.L);

end

Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt; 
Pd_MDL(jj)=length(find(Ns_MDL==num))./Nt;
Pd_AIC(jj)=length(find(Ns_AIC==num))./Nt;
Pd_IBIC(jj)=length(find(Ns_IBIC==num))./Nt;
Pd_ISSM(jj)=length(find(Ns_ISSM==num))./Nt;
Pd_MSRSE(jj)=length(find(Ns_MSRSE==num))./Nt;

end
%%
rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

xx=sensor_min:1:sensor_max;
hold on;
plot(xx,Pd_AIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(xx,Pd_MDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(xx,Pd_IBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(xx,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(xx,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');
plot(xx,Pd_MSRSE,'Color',rgbTriplet(6,:),'Marker','s');

xlabel('阵元数');
ylabel('正确检测概率');
axis([sensor_min sensor_max 0 1]);
legend('AIC','MDL','NBIC','GDE','ISSM','本文算法','Location','southeast');
% toc;