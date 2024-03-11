%% 将协方差矩阵矢量化，即用克罗内克积扩展阵列孔径，用空间平滑扩展阵列
% 参考论文：An Underdetermined Source Number Estimation Method for Non-Circular
% Targets Based on Sparse Array （2019）

function [Rv,subarray,z]=kronecker_spatial_smoothing(R)

[M,N] = size(R);
L = M*N;
z = zeros(1,L);

% 矢量化
for i=1:M
    z(1+(i-1)*N:i*N) = R(i,1:N);
end

% 空间平滑
subarray_length = L/2+1;
smooth_step = 1;
subarray_num = floor((L-subarray_length)/smooth_step)+1;
subarray = zeros(subarray_num,subarray_length);
for j=1:subarray_num
    subarray(j,:) = z(1+(j-1)*smooth_step:subarray_length+(j-1)*smooth_step);
end

Rv = subarray*subarray'/subarray_length;

% % 取平滑子阵列所有协方差的统计平均值
% Rv = zeros(subarray_length,subarray_length);
% for k=1:subarray_num
%     Rv = Rv + subarray(k,:)*subarray(k,:)'/subarray_length;
% end
% Rv = Rv/subarray_num;

end