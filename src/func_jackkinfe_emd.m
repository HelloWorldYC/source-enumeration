function [Ns_MIC,Ns_MSTDC] = func_jackkinfe_emd(signal,rou,times)

[signal_jackknife,index,L_jk] = Jackknife(signal,rou,times);
for i=1:times
    resignal = emd(signal_jackknife(i,:));
%   resignal = [signal_jackknife(i,:);resignal'];
    [M,~] = size(resignal);
    R = resignal*resignal'/L_jk;
    [~,v]=svd(R);
    T=diag(v);
    T1=T+sqrt(sum(T));
    [MIC,frequency_MIC(i)] = func_MIC(resignal,M,L_jk);
    [MSTDC,frequency_MSTDC(i)] = func_MSTDC(resignal,M,L_jk);
end

Ns_MIC = mode(frequency_MIC);
Ns_MSTDC = mode(frequency_MSTDC);
end