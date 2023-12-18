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
phi=[U1 U2 U3 U4 U5];
for j=1:5
teta(1,j)=pinv(phi(:,j))*y;
Y(:,j)=phi(:,j)*teta(1,j);
SSE(j)=sum(abs(Y(:,j)-y));
end
Y_n_1=y-Y(:,3);
phi_1=[U1 U2 U4 U5];
for j=1:4
teta_1(1,j)=pinv(phi_1(:,j))*Y_n_1;
Y_1(:,j)=phi_1(:,j)*teta_1(1,j);
SSE_1(j)=sum(abs(Y_1(:,j)-Y_n_1));
end
Y_n_2=Y_n_1-Y_1(:,1);
phi_2=[ U2 U4 U5];
for j=1:3
teta_2(1,j)=pinv(phi_2(:,j))*Y_n_2;
Y_2(:,j)=phi_2(:,j)*teta_2(1,j);
SSE_2(j)=sum(abs(Y_2(:,j)-Y_n_2));
end
%% %%%%%%%% %%%% result
% figure;
% plot(y,'linewidth',2);legend('sys ouput');grid on
PHI_N=[U1 U3 U4];
Condition_Number=cond(PHI_N)
disp('Forward selection [U1 U3 U4]' )
teta=pinv(PHI_N)*y;
y_hat=PHI_N*teta;
SSE_selected=sum((y-y_hat).^2)
figure
subplot(2,1,1)
plot(y,'b--','linewidth',2.5);grid on
hold on
plot(y_hat,'r-','linewidth',2);legend('real output','Forward select output')
subplot(2,1,2)
% error 
plot(abs(y_hat-y),'r-','linewidth',2);legend('ERROR')


%% other 3 regeressors
PHI_N1=[U1 U2 U3];
disp('[u1 u2 u3]')
teta1=pinv(PHI_N1)*y;
y_hat1=PHI_N1*teta1;
SSE_TOTAL1=sum((y-y_hat1).^2)
figure
plot(y,'r','linewidth',2);
hold on
plot(y_hat1,'b','linewidth',2);legend('sys ouput','random selected regs')