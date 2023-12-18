%% EXER_5_b_armax

%% System Identification Exercise_ first Series

%% Part 2_ARMAX
clc;
clear;
close all
% Generating Noise nad Input Signals
u=binrand([1:350],10,25,1,'normal');% generating input with use of T_binrand
varN=.01;% noise variance(average)
n=random('normal',0,varN,350,1);% generating noise with normal distribution N~(0,varN)
r =iddata([],u); % creates an iddata object for time-series data for input
v =iddata([],n); % creates an iddata object for time-series data for noise
%modeling
A = [1  -1.8006  0.8187];
B = [0  0.01813];
C=[1 -1.8187  0.8187];
sys= idpoly(A,B,C,1,1,varN,0.01)% generate Linear polynomial input-output model
yy = sim(sys,[r v]);%Simulate dynamic system with noise(v) and input(r)
yk = [yy,r];%Time domain data set with 350 samples
%Estimate parameters of ARX or AR model using least squares na=2 and nb=2(
%without any noise)
nb=2;nc=1% fixed nb
i=0;% diffrent na
for na=1:1:10;
    i=i+1;
    model_armax = armax(yk,[na nb nc 0]);% rders = [na nb nk] C=D=F=1
    Yhat = sim(model_armax,[r]);%Simulate dynamic system
    SSE=0;%Calculating SS
    SSET(i)=sum((yy.y-Yhat.y).^2);
end
%% plotting and sys estimatin
figure
na=1:1:10;
bar(na,SSET,'b','linewidth',3.5)
legend('SSER-na=0:1:10')
grid
[SSEmin,nmin]=min(SSET);
na=nmin
model_armax_best = armax(yk,[na nb nc 0])% rders = [na nb nk] C=D=F=1
Yhat = sim(model_armax_best,[r]);%Simulate dynamic system
SSE=0;%Calculating SS
SSE=sum((yy.y-Yhat.y).^2)
figure%ploting r 
subplot(2,2,1)
plot(r,'b','linewidth',2.5)
legend('INPUT')
grid on
%ploting real and model output
subplot(2,2,2)
plot(yy,'b','linewidth',90,Yhat,'r--','linewidth',2.5),
legend('real output','ARMAX model output')
grid on
%Specify discrete transfer functions 
discreteTF0=filt(sys.B,sys.F)*filt(1,sys.A);
discreteTFARMAX=filt(model_armax_best.B,model_armax_best.F)*filt(1,model_armax_best.A);
subplot(2,2,3)
pzplot(discreteTF0,'b')%Plot pole-zero map
legend('System')
subplot(2,2,4)
pzplot(discreteTFARMAX,'r')%Plot pole-zero map
legend('ARMAX model')