%% excel到format格式转换
clear;
clc;

%% 读入excel数据
name1 = 'E01_080_f14_waveform';
x1 = xlsread(['F:\信源数估计\code\MYC\main\labdata\',name1],'D100:G51000');

%% 输出format格式
