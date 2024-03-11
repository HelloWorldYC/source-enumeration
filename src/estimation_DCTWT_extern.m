%% ����˵��
%  ���ܣ� ����˫����С���ֽ��ⲿ��ͨ���źţ������ô�ͳ�㷨������Դ
%  ��Գ�� ��Ҷ��
%  �汾�� v1.0 - 2021.12.16

%% ���ⲿ���뵥ͨ���ź�
% ��ͨ���źŹ���Ƚϼ򵥣�����Ҫͨ������ʸ�����������죬��Ϊ��������Ԫ֮��Ĺ�ϵ
% ֱ�ӽ������ź�Դ�����������һ�𼴿ɣ��������������߶��ź�Դ�����棬�Ѿ�����
% �����䷽���Ӱ��
clear;
clc;
close all;
tic;

fs = 1024;
ts = 1/fs;
L = 1024;
t = (0:L-1)*ts;

% name = 'AD_Sample_nojam';
% signal = xlsread(['F:\��Դ������\code\��������\1801\',name],'A1:A1024');
% signal = importdata('3_80_fy_60_0du_4ch_cnr41_1.txt');
% signal = signal(1:L);

% 
% %% �Ե�ͨ���ź���DTCWT�ֽ�
% biort = 'near_sym_b';
% qshift = 'qshift_b';
% nlevel = 4;
% X = X';
% [Yl,Yh,Yscale] = dtwavexfm(X,nlevel,biort,qshift);
% z1 = dtwaveifm(Yl*0,Yh,biort,qshift,[1 0 0 0]);
% z2 = dtwaveifm(Yl*0,Yh,biort,qshift,[0 1 0 0]);
% z3 = dtwaveifm(Yl*0,Yh,biort,qshift,[0 0 1 0]);
% z4 = dtwaveifm(Yl*0,Yh,biort,qshift,[0 0 0 1]);
% z5 = dtwaveifm(Yl,Yh,biort,qshift,[0 0 0 0]);
% z = z1 + z2 + z3 + z4 + z5;
% 
% %%
% figure(2);
% subplot(5,1,1);
% plot(z1);
% subplot(5,1,2);
% plot(z2);
% subplot(5,1,3);
% plot(z3);
% subplot(5,1,4);
% plot(z4);
% subplot(5,1,5);
% plot(z5);
% figure(3);
% plot(z);
% 
% %% ���ɶ�ͨ���ź�
% % Y = [X z1 z2 z3 z4 z5];
% Y = [X z5 z4 z3 z2 z1];
% M = nlevel+2;
% Y = Y';
% R = Y*Y'/L;
% [GDE]=func_GDE(M,L,R);
% for k=1:M-1
%     if GDE(k)<0   
%         Ns_GDE=k-1; break;
%     end
% end


%% ��ͬ�ֽ�����µ�ʵ��
biort = 'near_sym_b';
qshift = 'qshift_b';
nlevel = 4:1:8;
accuracy_SNR = zeros(length(nlevel),1);
for level = 1:1:length(nlevel)
    unit = eye(nlevel(level));
    z = zeros(L,nlevel(level)+1);
    accuracy = 0;
    Ns_GDE = 0;
    X = signal;
%     X = X';
    [Yl,Yh,Yscale] = dtwavexfm(X,nlevel(level),biort,qshift);
    for j = 1:nlevel(level)
        z(:,j) = dtwaveifm(Yl*0,Yh,biort,qshift,unit(j,:));
    end
    z(:,nlevel(level)+1) = dtwaveifm(Yl,Yh,biort,qshift,zeros(1,nlevel(level)));
%     z = fliplr(z);
%     Y = [X z];
    Y = [z X];
    M = nlevel(level)+2;
    Y = Y';
    R = Y*Y'/L;
    [GDE]=func_GDE(M,L,R);
    for k=1:M-1
        if GDE(k)<0   
            Ns_GDE=k-1; break;
        end
    end
    accuracy_SNR(level) = Ns_GDE;
end
%% 
figure(3);
plot(nlevel,accuracy_SNR,'b^-');
xlabel('nlevel');
ylabel('estimation_number','Interpreter', 'none');
% title(strrep(name,'_','\_'));
title(name,'Interpreter', 'none');

toc;