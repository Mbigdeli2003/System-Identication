%% EXER_2_a_ARMAX
%%system identification exercises-first series
%generate clean input-output pairs, using system_2 the following input
%vector,for tmin=40 u=binrand([1:600];10;tmin;'normal')
%Try to fit an ARX model na=2,nb=2,nc=1 nk=0
clc
clear all 
u=binrand([1:600],10,40,1,'normal');


%A(q) = 1 - 1.72q^-1 + 0.9 q^-2 and B(q) = 0.48q^-1-0.48q^-2;C=D=f=1
%and Noise variance =0

A=[1 -1.72 0.9];
B=[0 0.48 -0.48];
sys= idpoly(A,B,1,1,1,0,0.1);
%create an iddata object
r=iddata([],u);
%simulate sys
yk = sim(sys,r);
yy= [yk,r];
%Estimate parameters of ARX or AR model using least squares na=2 and nb=2(
%without any noise)
model_armax=armax(yy,[2 2 1 0])%orders = [na nb nc nk] D=F=1
%Simulate dynamic system
Y_hat = sim(model_armax,r);
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
%ploting real and model output
figure 
plot(yk,'b','linewidth',3,Y_hat,'r-*')
legend('real output','model output ARMAX')
discreteTF0=filt(sys.B,sys.F)*filt(1,sys.A);
discreteTFARMAX=filt(model_armax.B,model_armax.F)*filt(1,model_armax.A);
figure;subplot(2,2,3)
pzplot(discreteTF0,'b')%Plot pole-zero map
legend('System')
subplot(2,2,4)
pzplot(discreteTFARMAX,'r')%Plot pole-zero map
legend('ARMAX model')

