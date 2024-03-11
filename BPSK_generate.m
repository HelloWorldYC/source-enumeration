%����˵��
%Rb����Ԫ����
%fc���ز�Ƶ��
%fs: ��������
%k�� ��Ԫ����
%A�� ��ֵ
function [code,bpsk] = BPSK_generate(Rb,fc,fs,k,A)
code = idinput(k,'rbs');  % ����α���������
code = code';
% code = ones(1,k);
N = k/Rb*fs;        % �ܹ��ж��ٸ�������
Npc = 1/Rb*fs;      % ÿ����Ԫ�Ĳ��������
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