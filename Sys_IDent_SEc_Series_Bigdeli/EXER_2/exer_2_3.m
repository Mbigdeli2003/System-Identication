clc
clear all
close all
n=1000;
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
phi=[ U1 U2 U3 U4 U5];
teta=pinv(phi)*y;
Y=phi*teta;
SSE=sum((Y-y).^2)
phi_1=[U1 U3 U4 U5];
teta_1=pinv(phi_1)*Y;
Y_1=phi_1*teta_1;
SSE_1=sum((y-Y_1).^2)
phi_2=[U1 U3 U4];
teta_2=pinv(phi_2)*Y_1;
Y_2=phi_2*teta_2;
SSE_2=sum((y-Y_2).^2)


%% %%%%%%%% %%%% result
% figure;
% plot(y,'linewidth',2);legend('sys ouput');grid on
PHI_N=[U1 U3 U5];
Condition_Number=cond(PHI_N)
disp('Forward selection [U1 U3 U4]' )
teta=pinv(PHI_N)*y;
y_hat=PHI_N*teta;
SSE_selected=sum((y-y_hat).^2)
figure
subplot(2,1,1)
plot(y,'b--','linewidth',2.5);grid on
hold on
plot(y_hat,'r-','linewidth',2);legend('real output','Backward Elimination output')
subplot(2,1,2)
% error 
plot(abs(y_hat-y),'r-','linewidth',2);legend('ERROR')







