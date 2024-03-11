%% 间隔抽样形成多通道
% x：输入要间隔抽样信号
% M: 要形成多通道的数目
function [Y]=interval_sampling(x,M)

N = length(x);
L = floor(N/M);
Y = zeros(M,L);

for i=1:M
    Y(i,:) = x(i:M:(M*L));
end

end