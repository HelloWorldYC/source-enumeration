%% 利用�?疏矩阵生成字典，再估�?
clear;
clc;
tic;

f0 = 15.48e4;
fs = 62e4;
fa = 2.3e4;
fb = 2.2e4;
Ns = 256;
L=Ns;
t=1:L;
num = 3; %信源�?

Array_Num=8;% 阵元�?
d=0.5; %线阵半径
lamda=1; %波长
kk=6;    %线阵
M=Array_Num;
num_max = 4;
% 入射�?
theta_jam=10:20:num_max*20;
degrad=pi/180;
%方位�?
alfa_jam=10:20:num_max*20;

for i=1:num_max
s_jam(i,:)=array_form(Array_Num,d,lamda,theta_jam(i),alfa_jam(i),kk);
end
% A=[s_jam(1,:);s_jam(2,:);s_jam(3,:)];%方向矩阵�?
A=[s_jam(1:num,:)];%方向矩阵�?
% A=[s_jam1];
A=A';

% �?疏表示参�?
param.L = 5;
param.K = 40;
param.numIteration = 50;
param.errorFlag = 0;
param.preserveDCAtom = 0;
% Dictionary = randn(L,param.K);
param.InitializationMethod = 'DataElements';
param.displayProgress = 0;

num_interval = L/num_max;
for i=1:num_max
[t1,at1,bt1,x(i,:)]=narrow_signal(fs,L,fa,fb,f0);
end
% x = randn(M,L);
% train_signal=x;
% x=[x1;x2;x3];
train_signal = zeros(M,L);
for j=1:num_max
temp=s_jam(1:j,:)'*x(1:j,:);
train_signal(:,(j-1)*num_interval+1:j*num_interval) = temp(:,1:num_interval);
end
rowrank = randperm(L);
train_signal = train_signal(:,rowrank);
[Dictionary_base,output_base] = KSVD(train_signal,param);
%%
Nt=50; %Monte次数
jj=0;
snr_min = -20;
snr_max = 20;
snr_length = snr_max-snr_min+1;

Pd_GDE=zeros(1,snr_length);
Pd_NBIC=zeros(1,snr_length);
Pd_ISSM=zeros(1,snr_length);
Pd_GDE_sparse=zeros(1,snr_length);
Pd_NBIC_sparse=zeros(1,snr_length);
Pd_ISSM_sparse=zeros(1,snr_length);
Pd_GDE_difference=zeros(1,snr_length);
Pd_NBIC_difference=zeros(1,snr_length);
Pd_ISSM_difference=zeros(1,snr_length);
Pd_GDE_DU=zeros(1,snr_length);
Pd_NBIC_DU=zeros(1,snr_length);
Pd_ISSM_DU=zeros(1,snr_length);
for SNR=snr_min:snr_max 
    disp(['SNR is ',num2str(SNR)]);
    Am=10^(SNR/10);
    jj=jj+1;
    Ns_GDE=zeros(1,Nt);
    Ns_NBIC=zeros(1,Nt);
    Ns_ISSM=zeros(1,Nt);
    Ns_GDE_sparse=zeros(1,Nt);
    Ns_NBIC_sparse=zeros(1,Nt);
    Ns_ISSM_sparse=zeros(1,Nt);
    Ns_GDE_difference=zeros(1,Nt);
    Ns_NBIC_difference=zeros(1,Nt);
    Ns_ISSM_difference=zeros(1,Nt);
    Ns_GDE_DU=zeros(1,Nt);
    Ns_NBIC_DU=zeros(1,Nt);
    Ns_ISSM_DU=zeros(1,Nt);

parfor cc=1:Nt
    x1 = zeros(num,L);
    for i=1:num
        [t1,at1,bt1,x1(i,:)]=narrow_signal(fs,L,fa,fb,f0);
    end
    signal=Am*x1;
    A1=A*signal; 
    X=awgn(A1,SNR,'measured');
    R=X*X'/L; %信号协方�?
    
    [GDE,Ns_GDE(cc)] = func_GDE(M,L,R);
    [IBIC,Ns_NBIC(cc)] = func_IBIC(1/(M*L),M,L,R);
    [ISSM,Ns_ISSM(cc)]=func_ISSM(X);
%     
%     [Dictionary,output] = KSVD(X,param);
    coef = OMP(Dictionary_base,X,param.L);
%     coef = full(output.CoefMatrix);
    coef = full(coef);
    coef_mean = sum(abs(coef))./sum(abs(coef)~=0);
    coef_re = full(coef);
    coef_re(abs(coef_re)<coef_mean)=0;
    [Uc,Sc,Vc]=svd(full(coef*coef'/L));
    Puc = Uc*(Uc'*Uc)^(-1)*Uc;
    Sc1 = diag(Sc);
    Sc1 = sqrt(sum(Sc1))+Sc1;
    Sc1 = diag(Sc1);
     
    resignal = Dictionary_base*coef_re;
    R2 = resignal*resignal'/L;
    [GDE_sparse,Ns_GDE_sparse(cc)] = func_GDE(M,L,R2);
    [NBIC_sparse,Ns_NBIC_sparse(cc)] = func_NBIC(1/(M*L),M,L,R2);
    [ISSM_sparse,Ns_ISSM_sparse(cc)]=func_ISSM(resignal);

    difference = resignal-X;
    Rd = difference*difference'/L;
    [GDE_difference,as1] = func_GDE(M,L,Rd);
    Ns_GDE_difference(cc) = M-as1;
    [NBIC_difference,as2] = func_NBIC(1/(M*L),M,L,Rd);
    Ns_NBIC_difference(cc) = M-as2;
    [ISSM_difference,as3]=func_ISSM(difference);
    Ns_ISSM_difference(cc) = M-as3;
    
    DU = Dictionary_base*Sc*Puc';
    [M2,L2]=size(DU);
    Rdu = DU*DU'/L2;
    [Udu,Sdu] = svd(Rdu);
    [GDE_DU,Ns_GDE_DU(cc)] = func_GDE(M2,L2,Rdu);
    [NBIC_DU,Ns_NBIC_DU(cc)] = func_NBIC(1/(M2*L2),M2,L2,Rdu);
    [ISSM_DU,Ns_ISSM_DU(cc)]=func_ISSM(DU);

end

Pd_GDE(jj)=length(find(Ns_GDE==num))./Nt;
Pd_NBIC(jj)=length(find(Ns_NBIC==num))./Nt;
Pd_ISSM(jj)=length(find(Ns_ISSM==num))./Nt;

Pd_GDE_sparse(jj)=length(find(Ns_GDE_sparse==num))./Nt;
Pd_NBIC_sparse(jj)=length(find(Ns_NBIC_sparse==num))./Nt;
Pd_ISSM_sparse(jj)=length(find(Ns_ISSM_sparse==num))./Nt;

Pd_GDE_difference(jj)=length(find(Ns_GDE_difference==num))./Nt;
Pd_NBIC_difference(jj)=length(find(Ns_NBIC_difference==num))./Nt;
Pd_ISSM_difference(jj)=length(find(Ns_ISSM_difference==num))./Nt;

Pd_GDE_DU(jj)=length(find(Ns_GDE_DU==num))./Nt;
Pd_NBIC_DU(jj)=length(find(Ns_NBIC_DU==num))./Nt;
Pd_ISSM_DU(jj)=length(find(Ns_ISSM_DU==num))./Nt;
end
 %%
xx=snr_min:snr_max;
figure;
plot(xx,Pd_GDE,'>-',xx,Pd_NBIC,'rs-',xx,Pd_ISSM,'b*-');
title(['白噪声下',num2str(Array_Num),'线阵估计',num2str(num),'个信�?']);
xlabel('不同信噪比（dB�?');
ylabel('正确�?测概�?(%)');
axis([snr_min snr_max 0 1]);
legend('GDE','IBIC','ISSM');

figure;
plot(xx,Pd_GDE_sparse,'>-',xx,Pd_NBIC_sparse,'rs-',xx,Pd_ISSM_sparse,'b*-');
title(['白噪声下',num2str(Array_Num),'线阵估计',num2str(num),'个信�?']);
xlabel('不同信噪比（dB�?');
ylabel('正确�?测概�?(%)');
axis([snr_min snr_max 0 1]);
legend('GDEsparse','IBICsparse','ISSMsparse');

figure;
plot(xx,Pd_GDE_difference,'>-',xx,Pd_NBIC_difference,'rs-',xx,Pd_ISSM_difference,'b*-');
title(['白噪声下',num2str(Array_Num),'线阵估计',num2str(num),'个信�?']);
xlabel('不同信噪比（dB�?');
ylabel('正确�?测概�?(%)');
axis([snr_min snr_max 0 1]);
legend('GDEdiff','IBICdiff','ISSMdiff');

figure;
plot(xx,Pd_GDE_DU,'>-',xx,Pd_NBIC_DU,'rs-',xx,Pd_ISSM_DU,'b*-');
title(['白噪声下',num2str(Array_Num),'线阵估计',num2str(num),'个信�?']);
xlabel('不同信噪比（dB�?');
ylabel('正确�?测概�?(%)');
axis([snr_min snr_max 0 1]);
legend('GDEdu','IBICdu','ISSMdu');
toc;
