%% ����˵��
%  ���ܣ� ����˫����С���ֽⵥͨ���źţ������ô�ͳ�㷨������Դ
%  ��Գ�� ��Ҷ��
%  �汾�� v1.0 - 2022.03.24

%%

clear;
clc;
close all;
tic;

fs = 5000;
L = 1024;
num = 3;    % ��Դ��
monte=100;           %Monte-Carloģ��Ĵ���
rou = 0.7;           % �۵���
times = 50;
level = 1;
biort = 'near_sym_b';
qshift = 'qshift_d';
SNR = -5:1:10;
z = cell(1,6);
resignal = zeros(6,L);
accuracy_SNR_AIC = zeros(1,length(SNR));
accuracy_SNR_MDL = zeros(1,length(SNR));
accuracy_SNR_GDE = zeros(1,length(SNR));
accuracy_SNR_BIC = zeros(1,length(SNR));

for i =1:length(SNR)
    Am=10^(SNR(i)/10);
    unit = eye(6);
    Ns_AIC=zeros(1,monte);
    Ns_MDL=zeros(1,monte);
    Ns_GDE=zeros(1,monte);
    Ns_BIC=zeros(1,monte);
    for mk=1:monte
        [t1,at1,bt1,s1]=narrow_signal(fs,L,100,120,1500);
        [t2,at2,bt2,s2]=narrow_signal(fs,L,80,95,2000);
        [t3,at3,bt3,s3]=narrow_signal(fs,L,130,120,1800);
        s = s1+s2+s3;
        noise = randn(1,L);
        x = Am*s+noise;
        x = x';
        [X,index] = Jackknife(x,rou,times);
       %% DCTWT�ֽ�
        %����ά��˫����С���������ֽⷽ�� ��15����45����75��������ζ�ŷֽ����㣬
        % ����һά��˵����û����ô��ֽⷽ��ģ����Էֽ������ѡ��
        [Yl,Yh,Yscale] = dtwavexfm2(X,level,biort,qshift);
        for j = 1:6
            z{j} = dtwaveifm2(Yl*0,Yh,biort,qshift,unit(:,j));
            [resignal(j,:),frequency] = inverse_Jackkinfe(z{j},index,L);
        end
        Y1 = resignal;
        [Mv,Lv] = size(resignal);
        Rv = Y1*Y1'/L;
%         [Rv,subarray,~]=kronecker_spatial_smoothing(R1);
%         [Mv,Lv] = size(Rv);
        [u1,v1] = svd(Rv);
        T1 = diag(v1);
        Tav1 = sqrt(sum(T1));
        T1 = T1+Tav1;
        [AIC,Ns_AIC(mk)] = func_AIC(Mv,Lv,T1);
        [MDL,Ns_MDL(mk)] = func_MDL(Mv,Lv,T1);
        [GDE,Ns_GDE(mk)] = func_GDE(Mv,Lv,Rv);
        [BIC,Ns_BIC(mk)] = func_BIC(1/(Mv*Lv),Mv,Lv,Rv);
    end
    accuracy_SNR_GDE(i)=length(find(Ns_GDE(:)==num))./monte;
    accuracy_SNR_MDL(i)=length(find(Ns_MDL(:)==num))./monte;
    accuracy_SNR_AIC(i)=length(find(Ns_AIC(:)==num))./monte;
    accuracy_SNR_BIC(i)=length(find(Ns_BIC(:)==num))./monte;
end
%% 
figure(1);
plot(SNR,accuracy_SNR_GDE,'b^-');
hold on;
plot(SNR,accuracy_SNR_MDL,'gX-');
plot(SNR,accuracy_SNR_AIC,'o-');
plot(SNR,accuracy_SNR_BIC,'p-');
xlabel('SNR');
ylabel('accuracy');
title('2D-DTCWT�ֽ��źŹ��ƽ��');
legend({'DTCWT_GDE','DTCWT_MDL','DTCWT_AIC','DTCWT_BIC'},'Interpreter','none');

figure(2);
subplot(3,1,1);plot(1:L,resignal(1,:));xlabel('t');ylabel('Amplitude');title('�ع��ź�1ʱ����');
subplot(3,1,2);plot(1:L,resignal(2,:));xlabel('t');ylabel('Amplitude');title('�ع��ź�2ʱ����');
subplot(3,1,3);plot(1:L,resignal(3,:));xlabel('t');ylabel('Amplitude');title('�ع��ź�3ʱ����');
figure(3);
subplot(3,1,1);plot(1:L,resignal(4,:));xlabel('t');ylabel('Amplitude');title('�ع��ź�4ʱ����');
subplot(3,1,2);plot(1:L,resignal(5,:));xlabel('t');ylabel('Amplitude');title('�ع��ź�5ʱ����');
subplot(3,1,3);plot(1:L,resignal(6,:));xlabel('t');ylabel('Amplitude');title('�ع��ź�6ʱ����');

figure(4);   
for i =1:3
   subplot(3,1,i);
   plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(resignal(i,:)))));
   xlabel('frequency/Hz');ylabel('Amplitude');title(['�ع��ź�',num2str(i),'Ƶ��']);
end
figure(5);
for i =1:3
   subplot(3,1,i);
   plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(resignal(i+3,:)))));
   xlabel('frequency/Hz');ylabel('Amplitude');title(['�ع��ź�',num2str(i+3),'Ƶ��']);
end

toc;