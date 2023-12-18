%% EXER_1_b
%%system identification exercises-first series

%Having generated input-ouput pairs %using the system 1, try to fit an ARX
%model with the parameters na=2, nb=2 and nk=0 to the data,three times and
%each for one the following input vectors.keep the generated data uncontaminated
clc
clear all
close all
u=randn(350,1);

%A(q) = 1 - 1.801 q^-1 + 0.8187 q^-2 and B(q) = 0.01813 q^-1;C=D=F=1
%and Noise variance =0
A= [1  -1.8006  0.8187];
B = [0  0.01813];
sys= idpoly(A,B,1,1,1,0,0.01);
%create an iddata object
r=iddata([],u);
%simulate sys
yk = sim(sys,r);
%Time domain data set with 350 samples.
yy= [yk,r];
%Estimate parameters of ARX or AR model using least squares na=2 and nb=2(
%without any noise)
model = arx(yy,[2 2 0])
%Simulate dynamic system
Y_hat = sim(model,r);
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
figure ;
subplot(2,2,1)
plot(yk,'b','linewidth',3)
legend('real output');
subplot(2,2,2)
plot(Y_hat,'r--','linewidth',2.5)
legend('model output')






