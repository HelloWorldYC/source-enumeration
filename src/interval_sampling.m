%% ��������γɶ�ͨ��
% x������Ҫ��������ź�
% M: Ҫ�γɶ�ͨ������Ŀ
function [Y]=interval_sampling(x,M)

N = length(x);
L = floor(N/M);
Y = zeros(M,L);

for i=1:M
    Y(i,:) = x(i:M:(M*L));
end

end