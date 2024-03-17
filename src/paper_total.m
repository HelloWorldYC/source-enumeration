%% �ܸ�m�ļ�
clear;
clc;
% paper1
% paper2
% paper3
% paper4
% paper5
% paper6
% paper7
% paper8
% paper9
% paper10
% paper11
% paper12
% paper13
% paper14
% paper15
% paper16

%% ��ͼ paper1
clear;
clc;
savefilename = './detection_probability/paper1_whitenoise_snr-20to20_snapshot256_sources3_sensors8.mat';
snr_min = -20;
snr_max = 20;

Pd_RAIC=load(savefilename,'Pd_AIC');
Pd_RAIC=Pd_RAIC.Pd_AIC;
Pd_RMDL=load(savefilename,'Pd_MDL');
Pd_RMDL=Pd_RMDL.Pd_MDL;
Pd_RNBIC=load(savefilename,'Pd_NBIC');
Pd_RNBIC=Pd_RNBIC.Pd_NBIC;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

xx=snr_min:snr_max;

hold on;
plot(xx,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(xx,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(xx,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');

box on;
grid on;
xlabel('�����(dB)');
ylabel('��ȷ������');
axis([snr_min snr_max 0 1]);
legend('AIC','MDL','NBIC','Location','southeast');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/�����°�������ʵ�鲻ͬ�����.png', '-dpng', '-r600');


%% ��ͼ paper2
clear;
clc;
savefilename = './detection_probability/paper2_whitenoise_snr0_snapshot30to300_sources3_sensors8.mat';
L_circle_min = 30;
L_circle_max = 300;

Pd_RAIC=load(savefilename,'Pd_AIC');
Pd_RAIC=Pd_RAIC.Pd_AIC;
Pd_RMDL=load(savefilename,'Pd_MDL');
Pd_RMDL=Pd_RMDL.Pd_MDL;
Pd_RNBIC=load(savefilename,'Pd_NBIC');
Pd_RNBIC=Pd_RNBIC.Pd_NBIC;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

L_circle = L_circle_min:10:L_circle_max;

hold on;
plot(L_circle,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(L_circle,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(L_circle,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');

box on;
grid on;
xlabel('������');
ylabel('��ȷ������');
axis([min(L_circle) max(L_circle) 0 1]);
legend('AIC','MDL','NBIC','Location','southeast');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/�����°�������ʵ�鲻ͬ������.png', '-dpng', '-r600');


%% ��ͼ paper3
clear;
clc;
savefilename = './detection_probability/paper3_whitenoise_snr10_snapshot256_sources1to6_sensors8.mat';
num_circle = 1:1:6;

Pd_RAIC=load(savefilename,'Pd_AIC');
Pd_RAIC=Pd_RAIC.Pd_AIC;
Pd_RMDL=load(savefilename,'Pd_MDL');
Pd_RMDL=Pd_RMDL.Pd_MDL;
Pd_RNBIC=load(savefilename,'Pd_NBIC');
Pd_RNBIC=Pd_RNBIC.Pd_NBIC;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

hold on;
plot(num_circle,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(num_circle,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(num_circle,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');

box on;
grid on;
xlabel('��Դ��');
ylabel('��ȷ������');
axis([min(num_circle) max(num_circle) 0 1]);
legend('AIC','MDL','NBIC','Location','southwest');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/�����°�������ʵ�鲻ͬ��Դ��.png', '-dpng', '-r600');


%% ��ͼ paper4
clear;
clc;
savefilename = './detection_probability/paper4_whitenoise_snr-6_snapshot256_sources3_sensors8to30.mat';

Pd_RAIC=load(savefilename,'Pd_AIC');
Pd_RAIC=Pd_RAIC.Pd_AIC;
Pd_RMDL=load(savefilename,'Pd_MDL');
Pd_RMDL=Pd_RMDL.Pd_MDL;
Pd_RNBIC=load(savefilename,'Pd_NBIC');
Pd_RNBIC=Pd_RNBIC.Pd_NBIC;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);
sensor_min = 8;
sensor_max = 30;
xx=sensor_min:1:sensor_max;

hold on;
plot(xx,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(xx,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(xx,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');

box on;
grid on;
xlabel('��Ԫ��');
ylabel('��ȷ������');
axis([sensor_min sensor_max 0 1]);
legend('AIC','MDL','NBIC','Location','southeast');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/�����°�������ʵ�鲻ͬ��Ԫ��.png', '-dpng', '-r600');


%% ��ͼ paper5
clear;
clc;
savefilename = './detection_probability/paper5_colornoise_snr-20to20_snapshot256_sources3_sensors8.mat';

snr_min = -20;
snr_max = 20;

Pd_RAIC=load(savefilename,'Pd_RAIC');
Pd_RAIC=Pd_RAIC.Pd_RAIC;
Pd_RMDL=load(savefilename,'Pd_RMDL');
Pd_RMDL=Pd_RMDL.Pd_RMDL;
Pd_RNBIC=load(savefilename,'Pd_RNBIC');
Pd_RNBIC=Pd_RNBIC.Pd_RNBIC;
Pd_GDE=load(savefilename,'Pd_GDE');
Pd_GDE=Pd_GDE.Pd_GDE;
Pd_ISSM=load(savefilename,'Pd_ISSM');
Pd_ISSM=Pd_ISSM.Pd_ISSM;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

xx=snr_min:snr_max;

hold on;
plot(xx,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(xx,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(xx,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(xx,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(xx,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');

box on;
grid on;
xlabel('�����(dB)');
ylabel('��ȷ������');
axis([snr_min snr_max 0 1]);
legend('RAIC','RMDL','RNBIC','GDE','ISSM','Location','southeast');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/������ɫ������ʵ�鲻ͬ�����.png', '-dpng', '-r600');


%% ��ͼ paper6
clear;
clc;
savefilename = './detection_probability/paper6_colornoise_snr4_snapshot30to300_sources3_sensors8.mat';

L_circle_min = 30;
L_circle_max = 300;
L_circle = L_circle_min:10:L_circle_max;

Pd_RAIC=load(savefilename,'Pd_RAIC');
Pd_RAIC=Pd_RAIC.Pd_RAIC;
Pd_RMDL=load(savefilename,'Pd_RMDL');
Pd_RMDL=Pd_RMDL.Pd_RMDL;
Pd_RNBIC=load(savefilename,'Pd_RNBIC');
Pd_RNBIC=Pd_RNBIC.Pd_RNBIC;
Pd_GDE=load(savefilename,'Pd_GDE');
Pd_GDE=Pd_GDE.Pd_GDE;
Pd_ISSM=load(savefilename,'Pd_ISSM');
Pd_ISSM=Pd_ISSM.Pd_ISSM;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);
 
hold on;
plot(L_circle,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(L_circle,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(L_circle,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(L_circle,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(L_circle,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');

box on;
grid on;
xlabel('������');
ylabel('��ȷ������');
axis([min(L_circle) max(L_circle) 0 1]);
legend('RAIC','RMDL','RNBIC','GDE','ISSM','Location','southeast');


% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/������ɫ������ʵ�鲻ͬ������.png', '-dpng', '-r600');


%% ��ͼ paper7
clear;
clc;
savefilename = './detection_probability/paper7_colornoise_snr10_snapshot256_sources1to6_sensors8.mat';

num_circle = 1:1:6;

Pd_RAIC=load(savefilename,'Pd_RAIC');
Pd_RAIC=Pd_RAIC.Pd_RAIC;
Pd_RMDL=load(savefilename,'Pd_RMDL');
Pd_RMDL=Pd_RMDL.Pd_RMDL;
Pd_RNBIC=load(savefilename,'Pd_RNBIC');
Pd_RNBIC=Pd_RNBIC.Pd_RNBIC;
Pd_GDE=load(savefilename,'Pd_GDE');
Pd_GDE=Pd_GDE.Pd_GDE;
Pd_ISSM=load(savefilename,'Pd_ISSM');
Pd_ISSM=Pd_ISSM.Pd_ISSM;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

hold on;
plot(num_circle,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(num_circle,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(num_circle,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(num_circle,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(num_circle,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');

box on;
grid on;
xlabel('��Դ��');
ylabel('��ȷ������');
axis([min(num_circle) max(num_circle) 0 1]);
legend('RAIC','RMDL','RNBIC','GDE','ISSM','Location','southwest');


% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/������ɫ������ʵ�鲻ͬ��Դ��.png', '-dpng', '-r600');



%% ��ͼ paper8
clear;
clc;
savefilename = './detection_probability/paper8_colornoise_snr0_snapshot256_sources3_sensors8to30.mat';

sensor_min = 8;
sensor_max = 30;

Pd_RAIC=load(savefilename,'Pd_RAIC');
Pd_RAIC=Pd_RAIC.Pd_RAIC;
Pd_RMDL=load(savefilename,'Pd_RMDL');
Pd_RMDL=Pd_RMDL.Pd_RMDL;
Pd_RNBIC=load(savefilename,'Pd_RNBIC');
Pd_RNBIC=Pd_RNBIC.Pd_RNBIC;
Pd_GDE=load(savefilename,'Pd_GDE');
Pd_GDE=Pd_GDE.Pd_GDE;
Pd_ISSM=load(savefilename,'Pd_ISSM');
Pd_ISSM=Pd_ISSM.Pd_ISSM;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

xx=sensor_min:1:sensor_max;

hold on;
plot(xx,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(xx,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(xx,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(xx,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(xx,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');

box on;
grid on;
xlabel('��Ԫ��');
ylabel('��ȷ������');
axis([sensor_min sensor_max 0 1]);
legend('RAIC','RMDL','RNBIC','GDE','ISSM','Location','southeast');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/������ɫ������ʵ�鲻ͬ��Ԫ��.png', '-dpng', '-r600');


%% ��ͼ paper9
clear;
clc;
savefilename = './detection_probability/paper9_whitenoise_snr-20to20_snapshot256_sources4_sensors8.mat';

snr_min = -20;
snr_max = 20;

Pd_RAIC=load(savefilename,'Pd_AIC');
Pd_RAIC=Pd_RAIC.Pd_AIC;
Pd_RMDL=load(savefilename,'Pd_MDL');
Pd_RMDL=Pd_RMDL.Pd_MDL;
Pd_RNBIC=load(savefilename,'Pd_NBIC');
Pd_RNBIC=Pd_RNBIC.Pd_NBIC;
Pd_GDE=load(savefilename,'Pd_GDE');
Pd_GDE=Pd_GDE.Pd_GDE;
Pd_ISSM=load(savefilename,'Pd_ISSM');
Pd_ISSM=Pd_ISSM.Pd_ISSM;
Pd_MSRSE=load(savefilename,'Pd_MSRSE');
Pd_MSRSE=Pd_MSRSE.Pd_MSRSE;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

xx=snr_min:snr_max;

hold on;
plot(xx,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(xx,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(xx,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(xx,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(xx,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');
plot(xx,Pd_MSRSE,'Color',rgbTriplet(6,:),'Marker','s');

box on;
grid on;
xlabel('�����(dB)');
ylabel('��ȷ������');
axis([snr_min snr_max 0 1]);
legend('AIC','MDL','NBIC','GDE','ISSM','�����㷨','Location','southeast');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/�����°�������ʵ�鲻ͬ�����.png', '-dpng', '-r600');


%% ��ͼ paper10
clear;
clc;
savefilename = './detection_probability/paper10_whitenoise_snr0_snapshot30to300_sources4_sensors8.mat';

L_circle_min = 30;
L_circle_max = 300;
L_circle = L_circle_min:10:L_circle_max;

Pd_RAIC=load(savefilename,'Pd_AIC');
Pd_RAIC=Pd_RAIC.Pd_AIC;
Pd_RMDL=load(savefilename,'Pd_MDL');
Pd_RMDL=Pd_RMDL.Pd_MDL;
Pd_RNBIC=load(savefilename,'Pd_NBIC');
Pd_RNBIC=Pd_RNBIC.Pd_NBIC;
Pd_GDE=load(savefilename,'Pd_GDE');
Pd_GDE=Pd_GDE.Pd_GDE;
Pd_ISSM=load(savefilename,'Pd_ISSM');
Pd_ISSM=Pd_ISSM.Pd_ISSM;
Pd_MSRSE=load(savefilename,'Pd_MSRSE');
Pd_MSRSE=Pd_MSRSE.Pd_MSRSE;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);
 
hold on;
plot(L_circle,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(L_circle,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(L_circle,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(L_circle,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(L_circle,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');
plot(L_circle,Pd_MSRSE,'Color',rgbTriplet(6,:),'Marker','s');

box on;
grid on;
xlabel('������');
ylabel('��ȷ������');
axis([min(L_circle) max(L_circle) 0 1]);
legend('AIC','MDL','NBIC','GDE','ISSM','�����㷨','Location','southeast');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/�����°�������ʵ�鲻ͬ������.png', '-dpng', '-r600');


%% ��ͼ paper11
clear;
clc;
savefilename = './detection_probability/paper11_whitenoise_snr10_snapshot256_sources1to6_sensors8.mat';

num_circle = 1:1:6;

Pd_RAIC=load(savefilename,'Pd_AIC');
Pd_RAIC=Pd_RAIC.Pd_AIC;
Pd_RMDL=load(savefilename,'Pd_MDL');
Pd_RMDL=Pd_RMDL.Pd_MDL;
Pd_RNBIC=load(savefilename,'Pd_NBIC');
Pd_RNBIC=Pd_RNBIC.Pd_NBIC;
Pd_GDE=load(savefilename,'Pd_GDE');
Pd_GDE=Pd_GDE.Pd_GDE;
Pd_ISSM=load(savefilename,'Pd_ISSM');
Pd_ISSM=Pd_ISSM.Pd_ISSM;
Pd_MSRSE=load(savefilename,'Pd_MSRSE');
Pd_MSRSE=Pd_MSRSE.Pd_MSRSE;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

hold on;
plot(num_circle,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(num_circle,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(num_circle,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(num_circle,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(num_circle,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');
plot(num_circle,Pd_MSRSE,'Color',rgbTriplet(6,:),'Marker','s');

box on;
grid on;
xlabel('��Դ��');
ylabel('��ȷ������');
axis([min(num_circle) max(num_circle) 0 1]);
legend('AIC','MDL','NBIC','GDE','ISSM','�����㷨','Location','southwest');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/�����°�������ʵ�鲻ͬ��Դ��.png', '-dpng', '-r600');


%% ��ͼ paper12
clear;
clc;
savefilename = './detection_probability/paper12_whitenoise_snr-4_snapshot256_sources4_sensors8to30.mat';

sensor_min = 8;
sensor_max = 30;

Pd_RAIC=load(savefilename,'Pd_AIC');
Pd_RAIC=Pd_RAIC.Pd_AIC;
Pd_RMDL=load(savefilename,'Pd_MDL');
Pd_RMDL=Pd_RMDL.Pd_MDL;
Pd_RNBIC=load(savefilename,'Pd_NBIC');
Pd_RNBIC=Pd_RNBIC.Pd_NBIC;
Pd_GDE=load(savefilename,'Pd_GDE');
Pd_GDE=Pd_GDE.Pd_GDE;
Pd_ISSM=load(savefilename,'Pd_ISSM');
Pd_ISSM=Pd_ISSM.Pd_ISSM;
Pd_MSRSE=load(savefilename,'Pd_MSRSE');
Pd_MSRSE=Pd_MSRSE.Pd_MSRSE;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

xx=sensor_min:1:sensor_max;
hold on;
plot(xx,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(xx,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(xx,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(xx,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(xx,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');
plot(xx,Pd_MSRSE,'Color',rgbTriplet(6,:),'Marker','s');

box on;
grid on;
xlabel('��Ԫ��');
ylabel('��ȷ������');
axis([sensor_min sensor_max 0 1]);
legend('AIC','MDL','NBIC','GDE','ISSM','�����㷨','Location','southeast');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/�����°�������ʵ�鲻ͬ��Ԫ��.png', '-dpng', '-r600');


%% ��ͼ paper13
clear;
clc;
savefilename = './detection_probability/paper13_colornoise_snr-20to20_snapshot256_sources4_sensors8.mat';

snr_min = -20;
snr_max = 20;

Pd_RAIC=load(savefilename,'Pd_RAIC');
Pd_RAIC=Pd_RAIC.Pd_RAIC;
Pd_RMDL=load(savefilename,'Pd_RMDL');
Pd_RMDL=Pd_RMDL.Pd_RMDL;
Pd_RNBIC=load(savefilename,'Pd_RNBIC');
Pd_RNBIC=Pd_RNBIC.Pd_RNBIC;
Pd_GDE=load(savefilename,'Pd_GDE');
Pd_GDE=Pd_GDE.Pd_GDE;
Pd_ISSM=load(savefilename,'Pd_ISSM');
Pd_ISSM=Pd_ISSM.Pd_ISSM;
Pd_MSRSE=load(savefilename,'Pd_MSRSE');
Pd_MSRSE=Pd_MSRSE.Pd_MSRSE;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

xx=snr_min:snr_max;

hold on;
plot(xx,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(xx,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(xx,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(xx,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(xx,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');
plot(xx,Pd_MSRSE,'Color',rgbTriplet(6,:),'Marker','s');
 
box on;
grid on;
xlabel('�����(dB)');
ylabel('��ȷ������');
axis([snr_min snr_max 0 1]);
legend('RAIC','RMDL','RNBIC','GDE','ISSM','�����㷨','Location','southeast');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/������ɫ������ʵ�鲻ͬ�����.png', '-dpng', '-r600');


%% ��ͼ paper14
clear;
clc;
savefilename = './detection_probability/paper14_colornoise_snr5_snapshot30to300_sources4_sensors8.mat';

L_circle_min = 30;
L_circle_max = 300;
L_circle = L_circle_min:10:L_circle_max;

Pd_RAIC=load(savefilename,'Pd_RAIC');
Pd_RAIC=Pd_RAIC.Pd_RAIC;
Pd_RMDL=load(savefilename,'Pd_RMDL');
Pd_RMDL=Pd_RMDL.Pd_RMDL;
Pd_RNBIC=load(savefilename,'Pd_RNBIC');
Pd_RNBIC=Pd_RNBIC.Pd_RNBIC;
Pd_GDE=load(savefilename,'Pd_GDE');
Pd_GDE=Pd_GDE.Pd_GDE;
Pd_ISSM=load(savefilename,'Pd_ISSM');
Pd_ISSM=Pd_ISSM.Pd_ISSM;
Pd_MSRSE=load(savefilename,'Pd_MSRSE');
Pd_MSRSE=Pd_MSRSE.Pd_MSRSE;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);
 
hold on;
plot(L_circle,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(L_circle,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(L_circle,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(L_circle,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(L_circle,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');
plot(L_circle,Pd_MSRSE,'Color',rgbTriplet(6,:),'Marker','s');

box on;
grid on;
xlabel('������');
ylabel('��ȷ������');
axis([min(L_circle) max(L_circle) 0 1]);
legend('RAIC','RMDL','RNBIC','GDE','ISSM','�����㷨','Location','southeast');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/������ɫ������ʵ�鲻ͬ������.png', '-dpng', '-r600');


%% ��ͼ paper15
clear;
clc;
savefilename = './detection_probability/paper15_colornoise_snr10_snapshot256_sources1to6_sensors8.mat';

num_circle = 1:1:6;

Pd_RAIC=load(savefilename,'Pd_RAIC');
Pd_RAIC=Pd_RAIC.Pd_RAIC;
Pd_RMDL=load(savefilename,'Pd_RMDL');
Pd_RMDL=Pd_RMDL.Pd_RMDL;
Pd_RNBIC=load(savefilename,'Pd_RNBIC');
Pd_RNBIC=Pd_RNBIC.Pd_RNBIC;
Pd_GDE=load(savefilename,'Pd_GDE');
Pd_GDE=Pd_GDE.Pd_GDE;
Pd_ISSM=load(savefilename,'Pd_ISSM');
Pd_ISSM=Pd_ISSM.Pd_ISSM;
Pd_MSRSE=load(savefilename,'Pd_MSRSE');
Pd_MSRSE=Pd_MSRSE.Pd_MSRSE;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

hold on;
plot(num_circle,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(num_circle,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(num_circle,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(num_circle,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(num_circle,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');
plot(num_circle,Pd_MSRSE,'Color',rgbTriplet(6,:),'Marker','s');

box on;
grid on;
xlabel('��Դ��');
ylabel('��ȷ������');
axis([min(num_circle) max(num_circle) 0 1]);
legend('RAIC','RMDL','RNBIC','GDE','ISSM','�����㷨','Location','southwest');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/������ɫ������ʵ�鲻ͬ��Դ��.png', '-dpng', '-r600');


%% ��ͼ paper16
clear;
clc;
savefilename = './detection_probability/paper16_colornoise_snr-5_snapshot256_sources4_sensors8to30.mat';

sensor_min = 8;
sensor_max = 30;

Pd_RAIC=load(savefilename,'Pd_RAIC');
Pd_RAIC=Pd_RAIC.Pd_RAIC;
Pd_RMDL=load(savefilename,'Pd_RMDL');
Pd_RMDL=Pd_RMDL.Pd_RMDL;
Pd_RNBIC=load(savefilename,'Pd_RNBIC');
Pd_RNBIC=Pd_RNBIC.Pd_RNBIC;
Pd_GDE=load(savefilename,'Pd_GDE');
Pd_GDE=Pd_GDE.Pd_GDE;
Pd_ISSM=load(savefilename,'Pd_ISSM');
Pd_ISSM=Pd_ISSM.Pd_ISSM;
Pd_MSRSE=load(savefilename,'Pd_MSRSE');
Pd_MSRSE=Pd_MSRSE.Pd_MSRSE;

rgbTriplet = 0.01*round(100*[062 043 109;...
    240 100 073;...
    255 170 050;...
    000 070 222;...
    046 158 43;...
    189 030 030]/255);

xx=sensor_min:1:sensor_max;
hold on;
plot(xx,Pd_RAIC,'Color',rgbTriplet(1,:),'Marker','*');
plot(xx,Pd_RMDL,'Color',rgbTriplet(2,:),'Marker','p');
plot(xx,Pd_RNBIC,'Color',rgbTriplet(3,:),'Marker','o');
plot(xx,Pd_GDE,'Color',rgbTriplet(4,:),'Marker','^');
plot(xx,Pd_ISSM,'Color',rgbTriplet(5,:),'Marker','d');
plot(xx,Pd_MSRSE,'Color',rgbTriplet(6,:),'Marker','s');

box on;
grid on;
xlabel('��Ԫ��');
ylabel('��ȷ������');
axis([sensor_min sensor_max 0 1]);
legend('RAIC','RMDL','RNBIC','GDE','ISSM','�����㷨','Location','southeast');

% ����ͼ�β�ָ�� DPI Ϊ 600
print('F:/�о�������/��ҵ���/��ҵ����/����ͼƬ/������ɫ������ʵ�鲻ͬ��Ԫ��.png', '-dpng', '-r600');
% toc;