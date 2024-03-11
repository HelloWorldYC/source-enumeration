%参数说明
%Rb：码元速率
%fc：载波频率
%fs: 采样速率
%k： 码元个数
%A： 幅值
function [code,bpsk] = BPSK_generate(Rb,fc,fs,k,A)
code = idinput(k,'rbs');  % 生成伪随机序列码
code = code';
% code = ones(1,k);
N = k/Rb*fs;        % 总共有多少个采样点数
Npc = 1/Rb*fs;      % 每个码元的采样点个数
l = 0;
bpsk = zeros(1,N);
for i=1:k
   for j = l:l+Npc-1
       if code(1,i)==-1
         bpsk(1,j+1) = A*cos(2*pi*fc*j/fs);
       elseif code(1,i)==1
         bpsk(1,j+1) = A*cos(2*pi*fc*j/fs + pi);
       end
   end
   l = l+Npc;
end