%% 论文第三章白噪声下不同快拍数实验
clear;
clc;
% tic;

f0 = 15.48e4;
fs = 62e4;
fa = 2.3e4;
fb = 2.2e4;
num=3; %信源数

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
A=s_jam(1:num,:);%方向矩阵；
A=A';

%%
Wp=2*pi*30;
Ws=2*pi*40;
Rp=0.5;
Rs=40;
fs1=100;
W=2*pi*fs1;
[N1,Wn]=buttord(2*Wp/W,2*Ws/W,Rp,Rs);
[b,a]=butter(N1,Wn);

Nt=2000; %Monte次数
jj=0;
snr = 0;
Am=10^(snr/10);
L_circle_min = 30;
L_circle_max = 300;
L_circle = L_circle_min:10:L_circle_max;
L_length = length(L_circle);
Pd_AIC=zeros(1,L_length);
Pd_MDL=zeros(1,L_length);
Pd_NBIC=zeros(1,L_length);
Pd_GDE=zeros(1,L_length);
Pd_ISSM=zeros(1,L_length);

for L=L_circle
    disp(['L is ',num2str(L)]);
    jj=jj+1;
    Ns_AIC=zeros(1,Nt);
    Ns_MDL=zeros(1,Nt);
    Ns_NBIC=zeros(1,Nt);
    Ns_GDE=zeros(1,Nt);
    Ns_ISSM=zeros(1,Nt);

for cc=1:Nt
    x1 = zeros(num,L);
    for i=1:num
        [t1,at1,bt1,x1(i,:)]=narrow_signal(fs,L,fa,fb,f0);
    end
    signal=Am*x1;
    A1=A*signal; 
    X=awgn(A1,snr,'measured');

    R=X*X'/L; %信号协方差

    [u,v]=svd(R);
    T=diag(v);
    [AIC,Ns_AIC(cc)] = func_AIC(M,L,T);
    [MDL,Ns_MDL(cc)] = func_MDL(M,L,T);
    [NBIC,Ns_NBIC(cc)] = func_NBIC(1/(M*L),M,L,R);
    [GDE,Ns_GDE(cc)] = func_GDE(M,L,R);
    [ISSM,Ns_ISSM(cc)]=func_ISSM(X);

end

Pd_AIC(jj)=length(find(Ns_AIC==num))./Nt;
Pd_MDL(jj)=length(find(Ns_MDL==num))./Nt;
Pd_NBIC(jj)=length(find(Ns_NBIC==num))./Nt;
Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt;
Pd_ISSM(jj)=length(find(Ns_ISSM==num))./Nt;

end

%%
savefilename = strcat('./detection_probability/paper2_whitenoise_snr', num2str(snr), ...
    '_snapshot',num2str(L_circle_min),'to',num2str(L_circle_max),'_sources',num2str(num),'_sensors',num2str(Array_Num),'.mat');
save(savefilename,'Pd_AIC','Pd_MDL','Pd_NBIC','Pd_GDE','Pd_ISSM');

 %%
 rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

hold on;
plot(L_circle,Pd_AIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(L_circle,Pd_MDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(L_circle,Pd_NBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(L_circle,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(L_circle,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');

box on;
grid on;
xlabel('快拍数');
ylabel('正确检测概率');
axis([min(L_circle) max(L_circle) 0 1]);
legend('AIC','MDL','NBIC','GDE','ISSM','Location','southeast');

% 保存图形并指定 DPI 为 600
print('F:/研究生事项/毕业答辩/毕业论文/论文图片/第三章白噪声下实验不同快拍数.png', '-dpng', '-r600');
% toc;
