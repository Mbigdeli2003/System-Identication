%% EXER_3_b_ARMAX
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
A = [1  -1.8006  0.8187];
B = [0  0.01813];
C=[1 -1.8187  0.8187];
n=random('normal',0,varN,350,1);% generating noise with normal distribution N~(0,varN)
sys= idpoly(A,B,C,1,1,varN,0.01)% generate Linear polynomial input-output model
r =iddata([],u); % creates an iddata object for time-series data for input
v =iddata([],n); % creates an iddata object for time-series data for nois
yk = sim(sys,[r v]);%Simulate dynamic system with noise(v) and input(r)
yy = [yk,r];%Time domain data set with 350 samples
%Estimate parameters of ARX or AR model using least squares na=2 and nb=2(
%without any noise)
model_armax=armax(yy,[2 2 1 0])%orders = [na nb nc nk] D=F=1
Y_hat = sim(model_armax,[r]);%Simulate dynamic system
S_E=0;
SE=sum((yk.y-Y_hat.y).^2);
SE=SE
%% plotting
%ploting r 
figure
plot(r,'b','linewidth',2.5)
legend('INPUT')
grid on
%ploting real and model output
figure 
plot(yk,'b','linewidth',3,Y_hat,'r-*');
legend('real output','model output ARMAX');
%Specify discrete transfer functions 
discreteTF0=filt(sys.B,sys.F)*filt(1,sys.A);
discreteTFARMAX=filt(model_armax.B,model_armax.F)*filt(1,model_armax.A);
figure;subplot(2,2,3)
pzplot(discreteTF0,'b')%Plot pole-zero map
legend('System')
subplot(2,2,4)
pzplot(discreteTFARMAX,'r')%Plot pole-zero map
legend('ARMAX model')


