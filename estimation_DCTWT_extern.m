%% 程序说明
%  功能： 利用双树复小波分解外部单通道信号，再利用传统算法估计信源
%  码猿： 马叶椿
%  版本： v1.0 - 2021.12.16

%% 从外部输入单通道信号
% 单通道信号构造比较简单，不需要通过方向矢量矩阵来构造，因为不存在阵元之间的关系
% 直接将几个信号源按增益叠加在一起即可，这个增益就是天线对信号源的增益，已经包括
% 了入射方向的影响
clear;
clc;
close all;
tic;

fs = 1024;
ts = 1/fs;
L = 1024;
t = (0:L-1)*ts;

% name = 'AD_Sample_nojam';
% signal = xlsread(['F:\信源数估计\code\海格数据\1801\',name],'A1:A1024');
% signal = importdata('3_80_fy_60_0du_4ch_cnr41_1.txt');
% signal = signal(1:L);

% 
% %% 对单通道信号做DTCWT分解
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
% %% 构成多通道信号
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


%% 不同分解层数下的实验
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