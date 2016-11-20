% matched Filter Code;
clc
clear all
close all

x1=1:12;
l=length(x1);
x2=ones(1,l);
x3=-x2;
n=0:l-1;
x4=sin(2*pi*n/l);


figure;
subplot(411);
stem(x1);
subplot(412);
stem(x2);
subplot(413);
stem(x3);
subplot(414)
stem(x4)

P1=sum(x1.^2);
P2=sum(x2.^2);
P3=sum(x3.^2);
P4=sum(x4.^2);

%Normalizing Signals
x1=x1./sqrt(P1);
x2=x2./sqrt(P2);
x3=x3./sqrt(P3);
x4=x4./sqrt(P4);

SNR=input('Enter The desired value of SNR in dBs:');
N0=10^(-SNR/10);

%generating noise signal of desired SNR
noise=sqrt(N0/l).*randn(1,l);

%Received Signal

r1=x1+noise;
r2=x2+noise;
r3=x3+noise;
r4=x4+noise;

%Matched Filter
M1=filter(fliplr(x1),1,r1);
M2=filter(fliplr(x2),1,r2);
M3=filter(fliplr(x3),1,r3);
M4=filter(fliplr(x4),1,r4);
figure;subplot(411);
plot(M1);
subplot(412);
plot(M2);
subplot(413);
plot(M3);
subplot(414);
plot(M4)

%correlation
X1=xcorr(r1,x1);
X2=xcorr(r2,x2);
X3=xcorr(r3,x3);
X4=xcorr(r4,x4);


figure;subplot(411);
plot(X1);
subplot(412);
plot(X2);
subplot(413);
plot(X3);
subplot(414);
 plot(X4)
