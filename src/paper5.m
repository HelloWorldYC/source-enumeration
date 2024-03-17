%% ���ĵ�����ɫ�����²�ͬ�����ʵ��
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
num=3; %��Դ��

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
Pd_RAIC=zeros(1,snr_length);
Pd_RMDL=zeros(1,snr_length);
Pd_RNBIC=zeros(1,snr_length);
Pd_GDE=zeros(1,snr_length);
Pd_ISSM=zeros(1,snr_length);
coef = cell(1,num_max);
for SNR=snr_min:snr_max 
    disp(['SNR is ',num2str(SNR)]);
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
    noise=randn(M,L);
    color_noise=filter(b,a,noise);        %�˲�������˹ɫ����
    X=A1+color_noise;

    R=X*X'/L; %�ź�Э����

    [u,v]=svd(R);
    T=diag(v);
    T1=T+sqrt(sum(T));
    [RAIC,Ns_RAIC(cc)] = func_AIC(M,L,T1);
    [RMDL,Ns_RMDL(cc)] = func_MDL(M,L,T1);
    [GDE,Ns_GDE(cc)] = func_GDE(M,L,R);
    [RNBIC,Ns_RNBIC(cc)] = func_RNBIC(1/(M*L),M,L,R);
    [ISSM,Ns_ISSM(cc)]=func_ISSM(X);
end

Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt; 
Pd_RMDL(jj)=length(find(Ns_RMDL==num))./Nt;
Pd_RAIC(jj)=length(find(Ns_RAIC==num))./Nt;
Pd_RNBIC(jj)=length(find(Ns_RNBIC==num))./Nt;
Pd_ISSM(jj)=length(find(Ns_ISSM==num))./Nt;

end

%%
savefilename = strcat('./detection_probability/paper5_colornoise_snr', num2str(snr_min), 'to', num2str(snr_max),...
    '_snapshot',num2str(L),'_sources',num2str(num),'_sensors',num2str(Array_Num),'.mat');
save(savefilename,'Pd_RAIC','Pd_RMDL','Pd_RNBIC','Pd_GDE','Pd_ISSM');

%%
rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

xx=snr_min:snr_max;

hold on;
plot(xx,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(xx,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(xx,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(xx,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(xx,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');

box on;
grid on;
xlabel('�����(dB)');
ylabel('��ȷ������');
axis([snr_min snr_max 0 1]);
legend('RAIC','RMDL','RNBIC','GDE','ISSM','Location','southeast');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/������ɫ������ʵ�鲻ͬ�����.png', '-dpng', '-r600');
% toc;