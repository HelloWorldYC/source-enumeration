clear;
clc;
close all;

p=1000;
n=1:p+1;
w=-pi:2*pi/1000:pi;
R=100; 
C=0.001;
wn=1/2*pi*R*C;
[b,a]=butter(1,wn); 
g=randn(1,1001); 
y=filter(b,a,g); 
at=y.*cos(w.*n); 
bt=y.*sin(w.*n); 
ft=at-bt; 

ct=y.*cos(w.*n); 
dt=y.*sin(w.*n); 
gt=ct-dt;

subplot(4,1,1);
plot(ft);
subplot(4,1,2);
ksdensity(ft);
subplot(4,1,3);
plot(abs(fft(ft)));
subplot(4,1,4);
ksdensity(gt);