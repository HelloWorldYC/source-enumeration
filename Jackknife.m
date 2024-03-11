%% Jackknife刀切法
% 码猿: 马叶椿
% 日期: 2022-01-14

%%
function [Y,index,L] = Jackknife(signal,rou,times)
% signal    输入的信号
% rou       折刀率
% times     折刀的次数
% Y         折刀后的信号

N = length(signal); % 输入信号的采样数
L = floor(rou*N);          % 折刀后的样本采样数
Y = zeros(times,L);
t = (1:L)*N/L;
t = floor(t);
index =zeros(times,L);
for i = 1:times
    T = randperm(N);
%     index = T(t);
%     index = sort(index);
%     Y(i,:) = signal(index);
    index(i,:) = T(1:L);
    index(i,:) = sort(index(i,:));
    Y(i,:) = signal(index(i,:));
end
end