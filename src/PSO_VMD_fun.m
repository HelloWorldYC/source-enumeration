%% PSO_VMD的待优化函数,用相关系数或方差来作为评价指标吧
function err=PSO_VMD_fun(range)

global snr;
nn = size(range);
alpha_range = range(:,1);
K_range = range(:,2);
n=nn(1);

fs = 10000;
L = 1000;
tau = 0;
[t1,at1,bt1,s1]=narrow_signal(fs,L,50,55,2000);
[t2,at2,bt2,s2]=narrow_signal(fs,L,60,65,2500);
[t3,at3,bt3,s3]=narrow_signal(fs,L,70,75,3000);
s = s1+s2+s3;
X = awgn(s,snr,'measured');
err=zeros(1,n);

for i=1:n
    [Y, u_hat, omega] = VMD(X, round(alpha_range(i)), tau, round(K_range(i)), false, 1, 1e-7);
    y1 = sum(Y,1);
    err(i) = var(y1-X); %用方差作为适应度
end

end