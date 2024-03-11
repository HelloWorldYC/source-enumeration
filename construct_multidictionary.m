%% 构造多个稀疏字典
function [Dictionary_base] = construct_multidictionary(fs,L,fa,fb,f0,param,num_max,s_jam,b,a)

[~,sensor]=size(s_jam);
x = zeros(num_max,L);
for i=1:num_max
[~,~,~,x(i,:)]=narrow_signal(fs,L,fa,fb,f0);
end

train_signal = cell(1,num_max);
% Pu_train_signal = cell(1,num_max);
Dictionary_base = cell(1,num_max);
% eigenvalue_nonoise = cell(1,num_max);
% U_nonoise = cell(1,num_max);
% Pu_train_signal = cell(1,num_max);
for j=1:num_max
train_signal{j}=s_jam(1:j,:)'*x(1:j,:);
train_signal{j}=awgn(train_signal{j},15,'measured');
% SNR=18;
% Am=10^(SNR/10);
% train_signal{j}=Am*train_signal{j};
% noise=randn(sensor,L);
% color_noise=filter(b,a,noise);        %滤波产生高斯色噪声
% train_signal{j}=train_signal{j}+color_noise;

% [U_nonoise{j},~,~]=svd(train_signal{j}*train_signal{j}'/L);
% Pu_train_signal{j} = U_nonoise{j}*((U_nonoise{j}'*U_nonoise{j})^(-1))*U_nonoise{j}';
[Dictionary_base{j},~] = KSVD(train_signal{j},param);
end

end