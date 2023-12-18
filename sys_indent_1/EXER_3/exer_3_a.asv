%% EXER_3_a_ARX
%%system identification exercises-first series
%considering the system_1 incited with the input(a)in the phase 1
%,contaminate the system with gaussian noise of v, with the mean of zero and variance of 2
%Try to model the system, using the five models of ARX,ARMAX, OE and BJ
%in the presence of minor,average,and major noise.compar each of the models
%outputs,parameters, and positions of poles and zeros with what you have
%from the system 1.whicj one of the models shows the performance among all
%of others
clc
clear all 
close all
u=binrand([1:350],10,25,1,'normal');
varN=2;% noise variance
n=random('normal',0,varN,350,1);% generating noise with normal distribution N~(0,varN)
r =iddata([],u); % creates an iddata object for time-series data for input
v =iddata([],n); % creates an iddata object for time-series data for noise
%modeling noise and sys
A = [1  -1.8006  0.8187];
B = [0  0.01813];
C=[1 -1.8187  0.8187];
sys= idpoly(A,B,C,1,1,varN,0.01)% generate Linear polynomial input-output model
yk= sim(sys,[r v]);%Simulate dynamic system with noise(v) and input(r)
yy= [yk,r];
%model estimation
model_arx = arx(yy,[2 2 0])
%Simulate dynamic system
Y_hat = sim(model_arx,r);
%Calculating SSE
S_E=0;
SE=sum((yk.y-Y_hat.y).^2);
SE=SE

%% plotting 
%ploting r 
figure
plot(r,'b','linewidth',2.5)
legend('INPUT')
grid on
%% ARX model
%ploting real and model output
figure 
plot(yk,'b','linewidth',3,Y_hat,'r-*')
legend('real output','model output ARX')
%Specify discrete transfer functions 
discreteTF0=filt(sys.B,sys.F)*filt(1,sys.A);
discreteTFARX=filt(model_arx.B,model_arx.F)*filt(1,model_arx.A);
figure;subplot(2,2,3)
pzplot(discreteTF0,'b')%Plot pole-zero map
legend('System')
subplot(2,2,4)
pzplot(discreteTFARX,'r')%Plot pole-zero map
legend('ARX model')