% 功能：将一个序列重构为延迟矩阵,即加窗处理
% 参数：
% x：输入序列
% delay_point：延迟点数，窗移动的步进
% column_num ：列数，即每一行元素的个数，或者说窗长
function delay_matrix = generate_delay_matrix(x,delay_point,column_num)
N = length(x);
row_num = ceil((N-column_num)/delay_point)+1;
delay_matrix = zeros(row_num-1,column_num);
for i = 1:row_num-1
    delay_matrix(i,:) = x(delay_point*(i-1)+1:delay_point*(i-1)+column_num);
end
