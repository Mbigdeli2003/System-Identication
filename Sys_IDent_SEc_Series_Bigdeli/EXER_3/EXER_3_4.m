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
y=zeros(N,1);
for i=3:N
    y(i)=-A(2:end)*y(i-1:-1:i-2)+B*u(i:-1:i-1)+C*e(i:-1:i-2);
end
z=iddata(y,u);
m_2 = n4sid(z,2)
figure
compare(z,m_2)
title('n=2')
m_4 = n4sid(z,4)
figure
compare(z,m_4)
title('n=4')
m_7 = n4sid(z,7)%,'foc','sim','ssp','can','ts',0)
figure
compare(z,m_7)
title('n=7')
m_9 = n4sid(z,9)%,'foc','sim','ssp','can','ts',0)
figure
compare(z,m_9)
title('n=9')

