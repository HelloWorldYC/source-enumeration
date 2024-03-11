%% Jackknife���з�
% ��Գ: ��Ҷ��
% ����: 2022-01-14

%%
function [Y,index,L] = Jackknife(signal,rou,times)
% signal    ������ź�
% rou       �۵���
% times     �۵��Ĵ���
% Y         �۵�����ź�

N = length(signal); % �����źŵĲ�����
L = floor(rou*N);          % �۵��������������
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