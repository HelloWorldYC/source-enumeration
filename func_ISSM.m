%% ISSM
% 码猿: 马叶椿
% 论文：Detection of the Number of Signals in Uniform Arrays by Invariants-Signal-Subspace Matching（2022）
% 版本：v1.0-2022.06.08

%%
function [ISSM,Ns_ISSM]=func_ISSM(Y)

[M,L] = size(Y);
ISSM = zeros(1,M-2);
% Y1 = Y(1:M-1,:);
% Y2 = Y(2:M,:);
R = Y*Y'/L;

% 第一对子阵列
R1 = R(1:M-1,:);
R2 = R(2:M,:);
% R1 = Y1*Y1'/L;
% R2 = Y2*Y2'/L;
[U1,~,~] = svd(R1);
[U2,~,~] = svd(R2);
% ISSM(1) = 2-(2*(U1(:,1)'*U2(:,1))^2)/((norm(U1(:,1),2))^2*(norm(U2(:,1),2))^2);
% I = eye(M-1);

% 第二对子阵列
R3 = R(1:M-2,:);
R4 = R(3:M,:);
[U3,~,~] = svd(R3);
[U4,~,~] = svd(R4);

for k=1:M-2
    u1_temp = U1(:,1:k);
    Puk1 = u1_temp*(u1_temp'*u1_temp)^(-1)*u1_temp';
%     IP1 = I-Puk11;
    u2_temp = U2(:,1:k);
    Puk2 = u2_temp*(u2_temp'*u2_temp)^(-1)*u2_temp';
    u3_temp = U3(:,1:k);
    Puk3 = u3_temp*(u3_temp'*u3_temp)^(-1)*u3_temp';
    u4_temp = U4(:,1:k);
    Puk4 = u4_temp*(u4_temp'*u4_temp)^(-1)*u4_temp';
%     IP2 = I-Puk12;
%     uhat1 = IP1*U1(:,k+1)/norm((IP1*U1(:,k+1)),2);
%     uhat2 = IP2*U2(:,k+1)/norm((IP2*U2(:,k+1)),2);
%     ISSM(k+1) = ISSM(k)+2-2*(norm(Puk1*uhat2,2))^2-2*(norm(Puk2*uhat1,2))^2-2*(uhat1'*U2(:,k+1))^2;
%     trace1 = trace(Puk1*Puk2);
    ISSM(k) = (2*k-2*trace(Puk1*Puk2)) + (2*k-2*trace(Puk3*Puk4));
end
[ISSM_min,Ns_ISSM] = min(ISSM);
% Ns_ISSM = Ns_ISSM-1;

end