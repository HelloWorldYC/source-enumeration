%% 用多通道的实验室数据检测其是否正常

clear;
clc;

SNR=-20:2:50;
accuracy_GDE=zeros(1,length(SNR));
accuracy_BIC=zeros(1,length(SNR));
num = 2;
k=0;
for snr=SNR
k=k+1;
name1 = ['EW2_',num2str(snr),'_waveform'];      
x1 = xlsread(['E:\MYC\main\labdata\',name1],'D101:D50100');
x2 = xlsread(['E:\MYC\main\labdata\',name1],'E101:E50100');
x3 = xlsread(['E:\MYC\main\labdata\',name1],'F101:F50100');
x4 = xlsread(['E:\MYC\main\labdata\',name1],'G101:G50100');


%%
L_total = length(x1);
L = 1000;
% M = 4;
times = floor(L_total/L);
Ns_GDE = zeros(1,times);
Ns_BIC = zeros(1,times);
level = 8;
unit = eye(level);
for i = 1:times
    index = (i-1)*L+1:i*L;
    signal = x2(index);
    [Yl,Yh,Yscale] = dtwavexfm(signal,level,'near_sym_b','qshift_b');%Yh中从高频到低频
    for j = 1:level
        z(:,j) = dtwaveifm(Yl,Yh,'near_sym_b','qshift_b',unit(j,:));
    end
    Y = [z(:,1:level) signal];
    Y = Y';
    [~,L1] = size(Y);
    Rv = Y*Y'/L1; 
    [M,n]=size(Rv);
    [GDE,Ns_GDE(i)]=func_GDE(M,n,Rv);
    [BIC,Ns_BIC(i)]=func_BIC(1/(M*n),M,n,Rv);
end
accuracy_GDE(k)=length(find(Ns_GDE==num))./times;
accuracy_BIC(k)=length(find(Ns_BIC==num))./times;
end
%%
figure(1);
plot(SNR,accuracy_GDE,'b^-',SNR,accuracy_BIC,'rx-');
legend('GDE','BIC');
xlabel('SNR');ylabel('accuracy');
title(['信源个数为',num2str(num),'    分解尺度为',num2str(level)]);
