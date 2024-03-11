% ���ܣ���һ�������ع�Ϊ�ӳپ���,���Ӵ�����
% ������
% x����������
% delay_point���ӳٵ��������ƶ��Ĳ���
% column_num ����������ÿһ��Ԫ�صĸ���������˵����
function delay_matrix = generate_delay_matrix(x,delay_point,column_num)
N = length(x);
row_num = ceil((N-column_num)/delay_point)+1;
delay_matrix = zeros(row_num-1,column_num);
for i = 1:row_num-1
    delay_matrix(i,:) = x(delay_point*(i-1)+1:delay_point*(i-1)+column_num);
end
