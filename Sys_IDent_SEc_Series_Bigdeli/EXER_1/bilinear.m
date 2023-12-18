%% bilinear
clc
clear all
close all
[u y]=input_output;
init1=[2 1.7 5 1 ];
%init1=[4 2 7 .5 ];
e=1*rand(1,300);
n=150;
teta_1=init1(1);
teta_2=init1(2);
teta_3=init1(3);
teta_4=init1(4);
%% 
for i=1:n
    k1(i)=y(i)-(4*u(i)^2);
    k2(i)=y(i)-(7*u(i)^0.5);
    y_1(i)=log(k1(i));
    y_2(i)=log(k2(i));
    
end

PHI=zeros(n,2);
P=zeros(2,2,n);
P(:,:,2)=10e+8*eye(2,2);
K=zeros(2,n);
theta_hat_12=zeros(2,n+1);
theta_hat_34=zeros(2,n+1);
for i=2:n
    phi=[1  log(u(i))]';
    K(:,i)=P(:,:,i)*phi/(1+phi'*P(:,:,i)*phi);
    P(:,:,i+1)=(eye(2,2)-K(:,i)*phi')*P(:,:,i);
    theta_hat_12(:,i+1)=theta_hat_12(:,i)+K(:,i)*(y_1(i)-phi'*theta_hat_12(:,i));
    theta_hat_34(:,i+1)=theta_hat_34(:,i)+K(:,i)*(y_2(i)-phi'*theta_hat_34(:,i));
    %teta(:,i)=[exp(theta_hat_12(1,i);theta_hat_12(2,i); exp(theta_hat_34(1,i); theta_hat_34(2,i)];
end
theta_RLS = theta_hat_12(:,end)';
theta_RLS_1 = theta_hat_34(:,end)';
teta=[exp(theta_RLS_1(1)) theta_RLS_1(2)   exp(theta_RLS(1)) theta_RLS(2)];

teta1=teta(1); teta2=teta(2);  teta3= teta(3); teta4= teta(4);

   %% ****** validation ******
for i=1:150
    y_hat(i)=teta1*u(i)^teta2+teta3*u(i)*(0.091*e(i))^teta4+e(i);
end

%% ****** test ******
for j=151:200
    y_hat1(j)=teta1*u(j)^teta2+teta3*u(j)*(0.091*e(j))^teta4;
end

figure
subplot(2,1,1)
plot(u)
legend('Input');
subplot(2,1,2)
plot(y , 'r')
legend('Output');

figure
subplot(2,1,1)
plot(y,'linewidth',2.5)
hold on
plot(y_hat,'r','linewidth',2.5)
legend('real output','Y train')
subplot(2,1,2)
plot(y(151:200),'linewidth',2.5)
hold on
plot(y_hat1(151:200),'r','linewidth',2.5)
legend('real ouput','Y Valid')
    
    
    