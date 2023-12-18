%% EXER_5_c_ararx_nb

%% System Identification Exercise_ first Series

%% Part 2_ARARX
clc;
clear;
close all
% Generating Noise nad Input Signals
u=binrand([1:350],10,25,1,'normal');% generating input with use of binrand
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
nc=1; nf=1; % fixed values in ARARX model
na=2; nd=2; i=0;
for nb=1:1:10;
    i=i+1;
    model_ararx =  pem(yk,[na nb nc nd nf 0]);% orders = [na nb nc nd nf nk] C=F=1
    Yhat = sim(model_ararx,[r]);%Simulate dynamic system
    SSE=0;%Calculating SS
    SSET(i)=sum((yy.y-Yhat.y).^2);
end
%% plotting and sys estimatin
figure
nb=1:1:10;
bar(nb,SSET,'b','linewidth',3.5)
legend('SSER-nb=0:1:10')
grid
[SSEmin,nmin]=min(SSET);
nb=nmin
model_ararx_best =  pem(yk,[na nb nc nd nf 0])% orders = [na nb nc nd nf nk] C=F=1
Yhat = sim(model_ararx_best,[r]);%Simulate dynamic system
SSE=0;%Calculating SS
SSE=sum((yy.y-Yhat.y).^2)
figure%ploting r 
subplot(2,2,1)
plot(r,'b','linewidth',2.5)
legend('INPUT')
grid on
%ploting real and model output
subplot(2,2,2)
plot(yy,'b','linewidth',90,Yhat,'r--','linewidth',2.5)
legend('real output','ARARX model output')
grid on
%Specify discrete transfer functions 
discreteTF0=filt(sys.B,sys.F)*filt(1,sys.A);
discreteTFARARX=filt(model_ararx_best.B,model_ararx_best.F)*filt(1,model_ararx_best.A);
subplot(2,2,3)
pzplot(discreteTF0,'b')%Plot pole-zero map
legend('System')
subplot(2,2,4)
pzplot(discreteTFARARX,'r')%Plot pole-zero map
legend('ARARX model')