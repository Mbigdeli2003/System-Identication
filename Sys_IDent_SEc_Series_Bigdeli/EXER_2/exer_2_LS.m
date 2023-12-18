clc
clear all
close all
n=1000;
A=[7 0.1 4 11 0.03];
%% input generation
%%U_1 input white noise mean=0;var=1
U1=random('normal',0,1,n,1);
%% U_2 input prbs
U2= idinput(n,'prbs');
%% U_3 input combination of u1 and u2
N=random('normal',0,(0.001)^(1/2),n,1);
U3=(0.2*U1)+(3*U2)+N;
%% U_4 input white noise mean=0.5 var=0.001
U4=random('normal',0.5,0.001^(1/2),n,1);
%% U_5% input cte=3
U5=3*ones(n,1);
%% white noise  
e=random('normal',0,0.1^(1/2),n,1);
%% output generation
y=7*U1+0.1*U2+4*U3+11*U4+0.03*U5+e;
phi=[U1 U2 U3 U4 U5 ];
teta=pinv(phi)*y;
disp('condition number of phi')
Cond_num=cond(phi)
figure;
plot(U1,'b','linewidth',1);grid on;hold on
plot(U2,'y','linewidth',1);grid on;hold on
plot(U3,'g','linewidth',1);grid on;hold on
plot(U4,'r','linewidth',1);grid on;hold on;
plot(U5,'y','linewidth',1);grid on;legend('u1','u2','u3','u4','u5')
%% y estimation
y_hat=teta(1)*U1+teta(2)*U2+teta(3)*U3+teta(4)*U4+teta(5)*U5;
figure;
plot(y,'b');legend('system output');grid on;

figure;
subplot(2,1,1)
plot(teta,'b-*','linewidth',2);legend('teta hat');grid on
subplot(2,1,2)
plot(teta'-A,'r-*','linewidth',2);legend('error parameter');grid on

%% estimated y_hat plot
figure;
subplot(2,1,1)
plot(y_hat,'b--','linewidth',2);hold on
plot(y,'r--','linewidth',1);legend('y hat','y');grid on
subplot(2,1,2)
plot(abs(y-y_hat),'r','linewidth',2);legend('error');grid on


















