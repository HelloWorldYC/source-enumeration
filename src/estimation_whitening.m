%% ��ͨ��ɫ���������£�ͨ���׻��˲�������תΪ���������ٹ�����Դ

clear;
clc;

Ns = 1024;
L=Ns;
t=1:L;
num=3; %��Դ��

Array_Num=8;% ��Ԫ��
d=0.5; %����뾶
lamda=1; %����
kk=6;    %����
M=Array_Num;
% �����
theta_jam1=10;
theta_jam2=40;
theta_jam3=70;
theta_jam4=70;
theta_jam5=88;
degrad=pi/180;
%��λ��
alfa_jam1=10;
alfa_jam2=50;
alfa_jam3=90;
alfa_jam4=140;
alfa_jam5=190;

s_jam1=array_form(Array_Num,d,lamda,theta_jam1,alfa_jam1,kk);
s_jam2=array_form(Array_Num,d,lamda,theta_jam2,alfa_jam2,kk);
s_jam3=array_form(Array_Num,d,lamda,theta_jam3,alfa_jam3,kk);
s_jam4=array_form(Array_Num,d,lamda,theta_jam4,alfa_jam4,kk);
s_jam5=array_form(Array_Num,d,lamda,theta_jam5,alfa_jam5,kk);
A=[s_jam1;s_jam2;s_jam3];%�������
% A=[s_jam1];
A=A';

% �����ͨ�˲���
Wp=2*pi*30;
Ws=2*pi*40;
Rp=0.5;
Rs=40;
fs=100;
W=2*pi*fs;
[N1,Wn]=buttord(2*Wp/W,2*Ws/W,Rp,Rs);
[b,a]=butter(N1,Wn);

Nt=100; %Monte����
jj=0;
Pd_GDE=zeros(1,41);
Pd_AIC=zeros(1,41);
Pd_MDL=zeros(1,41);
Pd_BIC=zeros(1,41);
Pd_BIC_qita10=zeros(1,41);
Pd_BIC_qita100=zeros(1,41);
for SNR=-30:10 
    Am=10^(SNR/10);
    jj=jj+1;
    Ns_AIC=zeros(1,Nt);
    Ns_MDL=zeros(1,Nt);
    Ns_GDE=zeros(1,Nt);
    Ns_BIC=zeros(1,Nt);
    Ns_BIC_qita10=zeros(1,Nt);
    Ns_BIC_qita100=zeros(1,Nt);
for cc=1:Nt
    x=randn(num,L);
    signal=Am*x;
    A1=A*signal; 
    noise=randn(M,L);
    color_noise=filter(b,a,noise);        %�˲�������˹ɫ����
    X=A1+color_noise;

    R=X*X'/L; %Դ�ź�Э����
    %%%%%%%%%%%%%�Ķ�������׻��˲���%%%%%%%%%%%%%%%
    Sn=color_noise*color_noise'/L;%������Э�������
    W=Sn^(-1/2);%����׻��˲���
    Y=W*R*W';%�׻��˲��任���Э�������
    %         Y=W*Rs*W';%�׻��˲��任���Э�������
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
[u,v]=svd(Y);
T=diag(v);
% T1=T+sqrt(sum(T));
[AIC,Ns_AIC(cc)] = func_AIC(M,L,T);
[MDL,Ns_MDL(cc)] = func_MDL(M,L,T);
[GDE,Ns_GDE(cc)] = func_GDE(M,L,Y);
[BIC,Ns_BIC(cc)] = func_BIC(1/(M*L),M,L,Y);
[BIC_qita10,Ns_BIC_qita10(cc)] = func_BIC(10/(M*L),M,L,Y);
[BIC_qita100,Ns_BIC_qita100(cc)] = func_BIC(100/(M*L),M,L,R);
end

Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt;
Pd_MDL(jj)=length(find(Ns_MDL==num))./Nt;
Pd_AIC(jj)=length(find(Ns_AIC==num))./Nt;
Pd_BIC(jj)=length(find(Ns_BIC==num))./Nt;
Pd_BIC_qita10(jj)=length(find(Ns_BIC_qita10==num))./Nt;
Pd_BIC_qita100(jj)=length(find(Ns_BIC_qita100==num))./Nt;
end
 %%
figure(1);
xx=-30:10;
plot(xx,Pd_AIC,'>-',xx,Pd_MDL,'rs-',xx,Pd_GDE,'o-',xx,Pd_BIC,'b^-');
title(['ɫ������',num2str(Array_Num),'�������',num2str(num),'����Դ']);
xlabel('��ͬ����ȣ�dB��');
ylabel('��ȷ������(%)');
axis([-30 10 0 1]);
legend('AIC','MDL','GDE','BIC');

figure(2);
xx=-30:10;
plot(xx,Pd_BIC,'>-',xx,Pd_BIC_qita10,'rs-',xx,Pd_BIC_qita100,'o-');
title(['ɫ������',num2str(Array_Num),'�������',num2str(num),'����Դ']);
xlabel('��ͬ����ȣ�dB��');
ylabel('��ȷ������(%)');
axis([-30 10 0 1]);
legend('qita:1','qita:10','qita:100');
%%