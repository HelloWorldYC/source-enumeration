%% ���ĵ����°������²�ͬ�����ʵ��
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
num=4; %��Դ��

Array_Num=8;% ��Ԫ��
d=0.5; %����뾶
lamda=1; %����
kk=6;    %����
M=Array_Num;
num_max = 6;
% �����
theta_jam=10:15:num_max*20;
degrad=pi/180;
%��λ��
alfa_jam=10:20:num_max*20;

s_jam = zeros(num_max,M);
for i=1:num_max
s_jam(i,:)=array_form(Array_Num,d,lamda,theta_jam(i),alfa_jam(i),kk);
end
A=s_jam(1:num,:);%�������
A=A';

% ϡ���ʾ����
param.L = 3;
param.K = 45;
param.numIteration = 50;
param.errorFlag = 1;
param.errorGoal = 1e-6;
param.preserveDCAtom = 0;
% Dictionary = randn(L,param.K);
param.InitializationMethod = 'DataElements';
param.displayProgress = 0;

% [Dictionary_base] = construct_multidictionary(fs,L,fa,fb,f0,param,num_max,s_jam);
% filename = strcat('./dictionaries/white/dictionary_white_sensor_', num2str(Array_Num), '.mat');
filename = strcat('dictionary_white.mat');
Dictionary_base = load(filename);
Dictionary_base = Dictionary_base.Dictionary_base;
%%
% �����ͨ�˲���
Wp=2*pi*30;
Ws=2*pi*35;
Rp=0.5;
Rs=40;
fs1=100;
W=2*pi*fs1;
[N1,Wn]=buttord(2*Wp/W,2*Ws/W,Rp,Rs);
[b,a]=butter(N1,Wn);

Nt=2000; %Monte����
jj=0;
snr_min = -20;
snr_max = 20;
snr_length = snr_max-snr_min+1;
Pd_GDE=zeros(1,snr_length);
Pd_AIC=zeros(1,snr_length);
Pd_MDL=zeros(1,snr_length);
Pd_NBIC=zeros(1,snr_length);
Pd_ISSM=zeros(1,snr_length);
Pd_MSRSE=zeros(1,snr_length);
coef = cell(1,num_max);
for SNR=snr_min:snr_max 
    disp(['SNR is ',num2str(SNR)]);
    Am=10^(SNR/10);
    jj=jj+1;
    Ns_AIC=zeros(1,Nt);
    Ns_MDL=zeros(1,Nt);
    Ns_GDE=zeros(1,Nt);
    Ns_NBIC=zeros(1,Nt);
    Ns_ISSM=zeros(1,Nt);
    Ns_MSRSE=zeros(1,Nt);
for cc=1:Nt
%     x=randn(num,L);
    x1 = zeros(num,L);
    for i=1:num
        [t1,at1,bt1,x1(i,:)]=narrow_signal(fs,L,fa,fb,f0);
    end
    signal=Am*x1;
    A1=A*signal; 
    X=awgn(A1,SNR,'measured');

    R=X*X'/L; %�ź�Э����

    [u,v]=svd(R);
    T=diag(v);
    [AIC,Ns_AIC(cc)] = func_AIC(M,L,T);
    [MDL,Ns_MDL(cc)] = func_MDL(M,L,T);
    [GDE,Ns_GDE(cc)] = func_GDE(M,L,R);
    [NBIC,Ns_NBIC(cc)] = func_NBIC(1/(M*L),M,L,R);
    [ISSM,Ns_ISSM(cc)]=func_ISSM(X);
    [MSRSE,Ns_MSRSE(cc)] = func_MSRSE(L,Dictionary_base,num_max,X,param.L);
end

Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt; 
Pd_MDL(jj)=length(find(Ns_MDL==num))./Nt;
Pd_AIC(jj)=length(find(Ns_AIC==num))./Nt;
Pd_NBIC(jj)=length(find(Ns_NBIC==num))./Nt;
Pd_ISSM(jj)=length(find(Ns_ISSM==num))./Nt;
Pd_MSRSE(jj)=length(find(Ns_MSRSE==num))./Nt;

end
%%
savefilename = strcat('./detection_probability/paper9_whitenoise_snr', num2str(snr_min), 'to', num2str(snr_max),...
    '_snapshot',num2str(L),'_sources',num2str(num),'_sensors',num2str(Array_Num),'.mat');
save(savefilename,'Pd_AIC','Pd_MDL','Pd_NBIC','Pd_GDE','Pd_ISSM','Pd_MSRSE');

%%
rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

xx=snr_min:snr_max;

hold on;
plot(xx,Pd_AIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(xx,Pd_MDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(xx,Pd_NBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(xx,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(xx,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');
plot(xx,Pd_MSRSE,'Color',rgbTriplet(6,:),'Marker','s');

box on;
grid on;
xlabel('�����(dB)');
ylabel('��ȷ������');
axis([snr_min snr_max 0 1]);
legend('AIC','MDL','NBIC','GDE','ISSM','�����㷨','Location','southeast');


% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/�����°�������ʵ�鲻ͬ�����.png', '-dpng', '-r600');
% toc;