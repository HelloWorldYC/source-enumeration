%% 将二维变成一维
% input：输入二维图像
% index：二维图像在一维中的索引
%     L：原一维信号的长度

function [signal,frequency]=inverse_Jackkinfe(input,index,L)
[m,n]=size(input);
signal = zeros(1,L);
frequency = zeros(1,L); % 源信号中各个时刻被提取次数

for i=1:m
    signal(index(i,:)) = signal(index(i,:)) + input(i,:);
    frequency(index(i,:)) = frequency(index(i,:)) + ones(1,n);
end

for j=1:L
    signal(j) = signal(j)/frequency(j);
end

end