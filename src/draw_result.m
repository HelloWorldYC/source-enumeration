%% ˫����С������

DTCWT_num = 6;

figure(1); % �ֽ���������ع����ǰ������Ƶ�ع��ź�
for i = 1:3
   subplot(3,1,i);
   plot(1:L,z(:,i));
   xlabel('t');ylabel('Amplitude');title(['�ع��ź�ʱ���Σ���Ƶ',num2str(i)]);
end
%%
figure(2); % �ֽ���������ع���ĺ�������Ƶ�ع��źż���Ƶ�ع��ź�
for i = 1:DTCWT_num-3
   subplot(DTCWT_num-3+1,1,i);
   plot(1:L,z(:,i+3));
   xlabel('t');ylabel('Amplitude');title(['�ع��ź�ʱ���Σ���Ƶ',num2str(i+3)]);
end
subplot(DTCWT_num-3+1,1,DTCWT_num-3+1);
plot(1:L,z(:,end));
xlabel('t');ylabel('Amplitude');title('�ع��ź�ʱ���Σ���Ƶ');

re = sum(z,2);  % �����ع��źŵ����
figure(3);
subplot(3,1,1);
plot(1:L,X);
xlabel('t');ylabel('Amplitude');title('��ͨ���ź�ʱ����');
subplot(3,1,2);
plot(1:L,re);
xlabel('t');ylabel('Amplitude');title('�ع��ź�ʱ����');
subplot(3,1,3);
plot(1:L,re-X);     % ԭʼ�ź����ع��źŲ�ֵ
xlabel('t');ylabel('Amplitude');title('��ֵ');

figure(4);     % ǰ������Ƶ�ع��ź�Ƶ��
for i =1:3
   subplot(3,1,i);
   plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(z(:,i)))));
   xlabel('frequency/Hz');ylabel('Amplitude');title(['�ع��ź�Ƶ�ף���Ƶ',num2str(i)]);
end

figure(5);
for i =1:DTCWT_num-3  % ���漸����Ƶ�ع��źż���Ƶ�ع��ź�Ƶ��
   subplot(DTCWT_num-3+1,1,i);
   plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(z(:,i+3)))));
   xlabel('frequency/Hz');ylabel('Amplitude');title(['�ع��ź��ź�Ƶ�ף���Ƶ',num2str(i+3)]);
end
subplot(DTCWT_num-3+1,1,DTCWT_num-3+1);
plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(z(:,end)))));
xlabel('frequency/Hz');ylabel('Amplitude');title('�ع��ź�Ƶ�ף���Ƶ');

figure(6);  % ԭʼ�ź�����ȫ�ع��ź�Ƶ��
subplot(2,1,1)
plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(X))));
xlabel('frequency/Hz');ylabel('Amplitude');title('��ͨ���ź�Ƶ��');
subplot(2,1,2)
plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(re))));
xlabel('frequency/Hz');ylabel('Amplitude');title('�ع��ź�Ƶ��');

figure(7);  % ԭʼ��ͨ���ź�����ȫ�ع��źŵĹ������ܶ�
subplot(2,1,1);
X_corr = xcorr(X);
plot((0:length(X_corr)-1)/length(X_corr)*fs-fs/2,abs(fftshift(fft(X_corr))));
xlabel('frequency/Hz');ylabel('Amplitude');title('��ͨ���źŹ������ܶ�');
subplot(2,1,2);
re_corr = xcorr(re);
plot((0:length(re_corr)-1)/length(re_corr)*fs-fs/2,abs(fftshift(fft(re_corr))));
xlabel('frequency/Hz');ylabel('Amplitude');title('�ع��źŹ������ܶ�');

figure(8);     % ǰ������Ƶ�ع��źŹ������ܶ�
for i =1:3
   subplot(3,1,i);
   temp = xcorr(z(:,i));
   plot((0:length(temp)-1)/length(temp)*fs-fs/2,abs(fftshift(fft(temp))));
   xlabel('frequency/Hz');ylabel('Amplitude');title(['�ع��źŹ������ܶȣ���Ƶ',num2str(i)]);
end

figure(9);
for i =1:DTCWT_num-3  % ���漸����Ƶ�ع��źż���Ƶ�ع��źŹ������ܶ�
   subplot(DTCWT_num-3+1,1,i);
   temp = xcorr(z(:,i+3));
   plot((0:length(temp)-1)/length(temp)*fs-fs/2,abs(fftshift(fft(temp))));
   xlabel('frequency/Hz');ylabel('Amplitude');title(['�ع��źŹ������ܶȣ���Ƶ',num2str(i+3)]);
end
subplot(DTCWT_num-3+1,1,DTCWT_num-3+1);
temp = xcorr(z(:,end));
plot((0:length(temp)-1)/length(temp)*fs-fs/2,abs(fftshift(fft(temp))));
xlabel('frequency/Hz');ylabel('Amplitude');title('�ع��źŹ������ܶȣ���Ƶ');
%%
figure(10);
subplot(3,1,1)
plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(s1))));
xlabel('frequency/Hz');ylabel('Amplitude');title('Դ�ź�s1Ƶ��');
subplot(3,1,2)
plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(s2))));
xlabel('frequency/Hz');ylabel('Amplitude');title('Դ�ź�s2Ƶ��');
subplot(3,1,3)
plot((0:L-1)/L*fs-fs/2,abs(fftshift(fft(s3))));
xlabel('frequency/Hz');ylabel('Amplitude');title('Դ�ź�s3Ƶ��');

figure(11);
subplot(3,1,1);
s1_corr = xcorr(s1);
plot((0:length(s1_corr)-1)/length(s1_corr)*fs-fs/2,abs(fftshift(fft(s1_corr))));
xlabel('frequency/Hz');ylabel('Amplitude');title('Դ�ź�s1�������ܶ�');
subplot(3,1,2);
s2_corr = xcorr(s2);
plot((0:length(s2_corr)-1)/length(s2_corr)*fs-fs/2,abs(fftshift(fft(s2_corr))));
xlabel('frequency/Hz');ylabel('Amplitude');title('Դ�ź�s2�������ܶ�');
subplot(3,1,3);
s3_corr = xcorr(s3);
plot((0:length(s3_corr)-1)/length(s3_corr)*fs-fs/2,abs(fftshift(fft(s3_corr))));
xlabel('frequency/Hz');ylabel('Amplitude');title('Դ�ź�s3�������ܶ�');

figure(12);
subplot(3,1,1);
plot(1:L,s1);
xlabel('t');ylabel('Amplitude');title('Դ�ź�s1ʱ����');
subplot(3,1,2);
plot(1:L,s2);
xlabel('t');ylabel('Amplitude');title('Դ�ź�s2ʱ����');
subplot(3,1,3);
plot(1:L,s3);
xlabel('t');ylabel('Amplitude');title('Դ�ź�s3ʱ����');

%%
figure(13);
for i =1:3
   subplot(3,1,i);
   plot((0:L-1)/L*fs-fs/2,angle(fftshift(fft(z(:,i)))));
   xlabel('frequency/Hz');ylabel('phase');title(['�ع��ź�Ƶ�ף���Ƶ',num2str(i)]);
end

figure(14);
for i =1:DTCWT_num-3  % ���漸����Ƶ�ع��źż���Ƶ�ع��ź�Ƶ��
   subplot(DTCWT_num-3+1,1,i);
   plot((0:L-1)/L*fs-fs/2,angle(fftshift(fft(z(:,i+3)))));
   xlabel('frequency/Hz');ylabel('phase');title(['�ع��ź��ź�Ƶ�ף���Ƶ',num2str(i+3)]);
end
subplot(DTCWT_num-3+1,1,DTCWT_num-3+1);
plot((0:L-1)/L*fs-fs/2,angle(fftshift(fft(z(:,end)))));
xlabel('frequency/Hz');ylabel('phase');title('�ع��ź�Ƶ�ף���Ƶ');