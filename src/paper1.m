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
Nt=2000; %Monte����
jj=0;
snr_min = -20;
snr_max = 20;
snr_length = snr_max-snr_min+1;
Pd_AIC=zeros(1,snr_length);
Pd_MDL=zeros(1,snr_length);
Pd_NBIC=zeros(1,snr_length);
coef = cell(1,num_max);
for SNR=snr_min:snr_max 
    disp(['SNR is ',num2str(SNR)]);
    Am=10^(SNR/10);
    jj=jj+1;
    Ns_AIC=zeros(1,Nt);
    Ns_MDL=zeros(1,Nt);
    Ns_NBIC=zeros(1,Nt);
for cc=1:Nt
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
    [NBIC,Ns_NBIC(cc)] = func_NBIC(1/(M*L),M,L,R);
end

Pd_MDL(jj)=length(find(Ns_MDL==num))./Nt;
Pd_AIC(jj)=length(find(Ns_AIC==num))./Nt;
Pd_NBIC(jj)=length(find(Ns_NBIC==num))./Nt;

end

%%
savefilename = strcat('./detection_probability/paper1_whitenoise_snr', num2str(snr_min), 'to', num2str(snr_max),...
    '_snapshot',num2str(L),'_sources',num2str(num),'_sensors',num2str(Array_Num),'.mat');
save(savefilename,'Pd_AIC','Pd_MDL','Pd_NBIC');

% pdaic=load(savefilename,'Pd_AIC');
% pdaic=pdaic.Pd_AIC;
% pdmdl=load(savefilename,'Pd_MDL');
% pdmdl=pdmdl.Pd_MDL;
% pdnbic=load(savefilename,'Pd_NBIC');
% pdnbic=pdnbic.Pd_NBIC;

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

box on;
grid on;
xlabel('�����(dB)');
ylabel('��ȷ������');
axis([snr_min snr_max 0 1]);
legend('AIC','MDL','NBIC','Location','southeast');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/�����°�������ʵ�鲻ͬ�����.png', '-dpng', '-r600');
% toc;