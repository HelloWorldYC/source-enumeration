%% ����ά���һά
% input�������άͼ��
% index����άͼ����һά�е�����
%     L��ԭһά�źŵĳ���

function [signal,frequency]=inverse_Jackkinfe(input,index,L)
[m,n]=size(input);
signal = zeros(1,L);
frequency = zeros(1,L); % Դ�ź��и���ʱ�̱���ȡ����

for i=1:m
    signal(index(i,:)) = signal(index(i,:)) + input(i,:);
    frequency(index(i,:)) = frequency(index(i,:)) + ones(1,n);
end

for j=1:L
    signal(j) = signal(j)/frequency(j);
end

end