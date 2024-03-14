%% 论文第三章色噪声下不同阵元数实验
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
num=3; %信源数

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
SNR = 0;
sensor_length = sensor_max-sensor_min + 1;
Pd_GDE=zeros(1,sensor_length);
Pd_RAIC=zeros(1,sensor_length);
Pd_RMDL=zeros(1,sensor_length);
Pd_RNBIC=zeros(1,sensor_length);
Pd_ISSM=zeros(1,sensor_length);

coef = cell(1,num_max);
for sensor=sensor_min:1:sensor_max

    disp(['ArrayNum is ',num2str(sensor)]);
    
    s_jam = zeros(num_max,sensor);
    for i=1:num_max
        s_jam(i,:)=array_form(sensor,d,lamda,theta_jam(i),alfa_jam(i),kk);
    end
    A=s_jam(1:num,:);%方向矩阵；
    
    A=A';
    
    Am=10^(SNR/10);
    jj=jj+1;
    Ns_RAIC=zeros(1,Nt);
    Ns_RMDL=zeros(1,Nt);
    Ns_GDE=zeros(1,Nt);
    Ns_RNBIC=zeros(1,Nt);
    Ns_ISSM=zeros(1,Nt);

for cc=1:Nt
    x1 = zeros(num,L);
    for i=1:num
        [t1,at1,bt1,x1(i,:)]=narrow_signal(fs,L,fa,fb,f0);
    end
    signal=Am*x1;
    A1=A*signal; 
    noise=randn(sensor,L);
    color_noise=filter(b,a,noise);        %滤波产生高斯色噪声
    X=A1+color_noise;

    R=X*X'/L; %信号协方差

    [u,v]=svd(R);
    T=diag(v);
    T1=T+sqrt(sum(T));
    [RAIC,Ns_RAIC(cc)] = func_AIC(sensor,L,T1);
    [RMDL,Ns_RMDL(cc)] = func_MDL(sensor,L,T1);
    [GDE,Ns_GDE(cc)] = func_GDE(sensor,L,R);
    [RNBIC,Ns_RNBIC(cc)] = func_RNBIC(1/(sensor*L),sensor,L,R);
    [ISSM,Ns_ISSM(cc)]=func_ISSM(X);

end

Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt; 
Pd_RMDL(jj)=length(find(Ns_RMDL==num))./Nt;
Pd_RAIC(jj)=length(find(Ns_RAIC==num))./Nt;
Pd_RNBIC(jj)=length(find(Ns_RNBIC==num))./Nt;
Pd_ISSM(jj)=length(find(Ns_ISSM==num))./Nt;

end

%%
savefilename = strcat('./detection_probability/paper8_colornoise_snr', num2str(SNR), ...
    '_snapshot',num2str(L),'_sources',num2str(num),'_sensors',num2str(sensor_min),'to',num2str(sensor_max),'.mat');
save(savefilename,'Pd_RAIC','Pd_RMDL','Pd_RNBIC','Pd_GDE','Pd_ISSM');

%%
rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

xx=sensor_min:1:sensor_max;

hold on;
plot(xx,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(xx,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(xx,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(xx,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(xx,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');

box on;
xlabel('阵元数');
ylabel('正确检测概率');
axis([sensor_min sensor_max 0 1]);
legend('RAIC','RMDL','RNBIC','GDE','ISSM','Location','southeast');


% 保存图形并指定 DPI 为 600
print('F:/研究生事项/毕业答辩/毕业论文/论文图片/第三章色噪声下实验不同阵元数.png', '-dpng', '-r600');
% toc;