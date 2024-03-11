%% 双树复小波部分

DTCWT_num = 6;

figure(1); % 分解出来分量重构后的前三个高频重构信号
for i = 1:3
   subplot(3,1,i);
   plot(1:L,z(:,i));
   xlabel('t');ylabel('Amplitude');title(['重构信号时域波形：高频',num2str(i)]);
end
%%
figure(2); % 分解出来分量重构后的后三个高频重构信号及低频重构信号
for i = 1:DTCWT_num-3
   subplot(DTCWT_num-3+1,1,i);
   plot(1:L,z(:,i+3));
   xlabel('t');ylabel('Amplitude');title(['重构信号时域波形：高频',num2str(i+3)]);
end
subplot(DTCWT_num-3+1,1,DTCWT_num-3+1);
plot(1:L,z(:,end));
xlabel('t');ylabel('Amplitude');title('重构信号时域波形：低频');

re = sum(z,2);  % 所有重构信号的相加
figure(3);
subplot(3,1,1);
plot(1:L,X);
xlabel('t');ylabel('Amplitude');title('单通道信号时域波形');
subplot(3,1,2);
plot(1:L,re);
xlabel('t');ylabel('Amplitude');title('重构信号时域波形');
subplot(3,1,3);
plot(1:L,re-X);     % 原始信号与重构信号差值
xlabel('t');ylabel('Amplitude');title('差值');

figure(4);     % 前三个高频重构信号频谱
for i =1:3
   subplot(3,1,i);
   plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(z(:,i)))));
   xlabel('frequency/Hz');ylabel('Amplitude');title(['重构信号频谱：高频',num2str(i)]);
end

figure(5);
for i =1:DTCWT_num-3  % 后面几个高频重构信号及低频重构信号频谱
   subplot(DTCWT_num-3+1,1,i);
   plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(z(:,i+3)))));
   xlabel('frequency/Hz');ylabel('Amplitude');title(['重构信号信号频谱：高频',num2str(i+3)]);
end
subplot(DTCWT_num-3+1,1,DTCWT_num-3+1);
plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(z(:,end)))));
xlabel('frequency/Hz');ylabel('Amplitude');title('重构信号频谱：低频');

figure(6);  % 原始信号与完全重构信号频谱
subplot(2,1,1)
plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(X))));
xlabel('frequency/Hz');ylabel('Amplitude');title('单通道信号频谱');
subplot(2,1,2)
plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(re))));
xlabel('frequency/Hz');ylabel('Amplitude');title('重构信号频谱');

figure(7);  % 原始单通道信号与完全重构信号的功率谱密度
subplot(2,1,1);
X_corr = xcorr(X);
plot((0:length(X_corr)-1)/length(X_corr)*fs-fs/2,abs(fftshift(fft(X_corr))));
xlabel('frequency/Hz');ylabel('Amplitude');title('单通道信号功率谱密度');
subplot(2,1,2);
re_corr = xcorr(re);
plot((0:length(re_corr)-1)/length(re_corr)*fs-fs/2,abs(fftshift(fft(re_corr))));
xlabel('frequency/Hz');ylabel('Amplitude');title('重构信号功率谱密度');

figure(8);     % 前三个高频重构信号功率谱密度
for i =1:3
   subplot(3,1,i);
   temp = xcorr(z(:,i));
   plot((0:length(temp)-1)/length(temp)*fs-fs/2,abs(fftshift(fft(temp))));
   xlabel('frequency/Hz');ylabel('Amplitude');title(['重构信号功率谱密度：高频',num2str(i)]);
end

figure(9);
for i =1:DTCWT_num-3  % 后面几个高频重构信号及低频重构信号功率谱密度
   subplot(DTCWT_num-3+1,1,i);
   temp = xcorr(z(:,i+3));
   plot((0:length(temp)-1)/length(temp)*fs-fs/2,abs(fftshift(fft(temp))));
   xlabel('frequency/Hz');ylabel('Amplitude');title(['重构信号功率谱密度：高频',num2str(i+3)]);
end
subplot(DTCWT_num-3+1,1,DTCWT_num-3+1);
temp = xcorr(z(:,end));
plot((0:length(temp)-1)/length(temp)*fs-fs/2,abs(fftshift(fft(temp))));
xlabel('frequency/Hz');ylabel('Amplitude');title('重构信号功率谱密度：低频');
%%
figure(10);
subplot(3,1,1)
plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(s1))));
xlabel('frequency/Hz');ylabel('Amplitude');title('源信号s1频谱');
subplot(3,1,2)
plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(s2))));
xlabel('frequency/Hz');ylabel('Amplitude');title('源信号s2频谱');
subplot(3,1,3)
plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(s3))));
xlabel('frequency/Hz');ylabel('Amplitude');title('源信号s3频谱');

figure(11);
subplot(3,1,1);
s1_corr = xcorr(s1);
plot((0:length(s1_corr)-1)/length(s1_corr)*fs-fs/2,abs(fftshift(fft(s1_corr))));
xlabel('frequency/Hz');ylabel('Amplitude');title('源信号s1功率谱密度');
subplot(3,1,2);
s2_corr = xcorr(s2);
plot((0:length(s2_corr)-1)/length(s2_corr)*fs-fs/2,abs(fftshift(fft(s2_corr))));
xlabel('frequency/Hz');ylabel('Amplitude');title('源信号s2功率谱密度');
subplot(3,1,3);
s3_corr = xcorr(s3);
plot((0:length(s3_corr)-1)/length(s3_corr)*fs-fs/2,abs(fftshift(fft(s3_corr))));
xlabel('frequency/Hz');ylabel('Amplitude');title('源信号s3功率谱密度');

figure(12);
subplot(3,1,1);
plot(1:L,s1);
xlabel('t');ylabel('Amplitude');title('源信号s1时域波形');
subplot(3,1,2);
plot(1:L,s2);
xlabel('t');ylabel('Amplitude');title('源信号s2时域波形');
subplot(3,1,3);
plot(1:L,s3);
xlabel('t');ylabel('Amplitude');title('源信号s3时域波形');

%%
figure(13);
for i =1:3
   subplot(3,1,i);
   plot((0:L-1)/L*fs-fs/2,angle(fftshift(fft(z(:,i)))));
   xlabel('frequency/Hz');ylabel('phase');title(['重构信号频谱：高频',num2str(i)]);
end

figure(14);
for i =1:DTCWT_num-3  % 后面几个高频重构信号及低频重构信号频谱
   subplot(DTCWT_num-3+1,1,i);
   plot((0:L-1)/L*fs-fs/2,angle(fftshift(fft(z(:,i+3)))));
   xlabel('frequency/Hz');ylabel('phase');title(['重构信号信号频谱：高频',num2str(i+3)]);
end
subplot(DTCWT_num-3+1,1,DTCWT_num-3+1);
plot((0:L-1)/L*fs-fs/2,angle(fftshift(fft(z(:,end)))));
xlabel('frequency/Hz');ylabel('phase');title('重构信号频谱：低频');