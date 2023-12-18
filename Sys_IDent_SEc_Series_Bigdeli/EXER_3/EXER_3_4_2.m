%% Exe 3 part 5
clc;    clear all;  close all

N=150;
Ts=0.01;
u= idinput(N,'prbs');
e=random('normal',0,0.2,N,1);
%% system parameters
A=[1 0.6 0.08];
B=[0.5 -0.9];
C=[3 6 3];
%% initialization
varN=2;% noise variance
A = [1  -1.8006  0.8187];
B = [0  0.01813];
C=[1 -1.8187  0.8187];
n=random('normal',0,varN,350,1);% generating noise with normal distribution N~(0,varN)
sys= idpoly(A,B,C,1,1,varN,0.01)% generate Linear polynomial input-output model
r =iddata([],u); % creates an iddata object for time-series data for input
v =iddata([],e); % creates an iddata object for time-series data for nois
yk = sim(sys,[r v]);%Simulate dynamic system with noise(v) and input(r)
yy = [yk,r];%Time domain data set with 350 samples
y=zeros(N,1);
for i=3:N
    y(i)=-A(2:end)*y(i-1:-1:i-2)+B*u(i:-1:i-1)+C*e(i:-1:i-2);
end
z=iddata(y,u);
m_2 = n4sid(z,2);
orders=[2 2];
m2_a =armax(yy,[2 2 1 0])%orders = [na nb nc nk] D=F=1
Y_R_2 = sim(m2_a,[r v])
ARMAX_2=iddata(Y_R_2.y,u);
figure
compare(m_2,ARMAX_2)
title('n=2')
m_4 = n4sid(z,4);
orders=[2 2];
m2_a =armax(yy,[4 4 2 0])%orders = [na nb nc nk] D=F=1
Y_R_2 = sim(m2_a,[r v])
ARMAX_4=iddata(Y_R_2.y,u);
figure
compare(m_4,ARMAX_4)
title('n=4')
m_7 = n4sid(z,7);
m2_a =armax(yy,[7 7 2 0])%orders = [na nb nc nk] D=F=1
Y_R_2 = sim(m2_a,[r v])
ARMAX_7=iddata(Y_R_2.y,u);
figure
compare(m_7,ARMAX_7)
title('n=7')
% m_9 = n4sid(z,9)%,'foc','sim','ssp','can','ts',0)
% m9_a = armax(z,9);
% figure
% compare(m9_a,m_9)
% title('n=9')
% 
