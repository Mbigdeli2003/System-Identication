clc
clear all
close all
%% input
%% number of datas
n=200
U= idinput(n,'prbs');
%% White noise generation
e=random('normal',0,0.2^(1/2),n,1);
%% system parameters
A=[1 0.6 0.08];
B=[0.5 -0.9];
C=[3 6 3];
y=zeros(n,1);
phi= zeros(1,6);
k = zeros(length(phi(1,:)));
p = zeros(length(phi(1,:)));
alpha =1e12;
p(:,:,2) = alpha*eye(length(phi(1,:)));
teta = [zeros(6,1),zeros(6,1),zeros(6,1)];
E=ones(1,3);
y_hat=zeros(1,n);
for i=3:150
    y(i)=-(A(2)*y(i-1)+A(3)*y(i-2))+B(1)*U(i)+B(2)*U(i-1)+C(1)*e(i)+C(2)*e(i-1)+C(3)*e(i-2);
      phi(i,:)=[-y(i-1) -y(i-2) U(i) U(i-1) e(i-1) e(i-2)];
    k(:,i)=p(:,:,i-1)*phi(i,:)'/(1+phi(i,:)*p(:,:,i-1)*phi(i,:)');
    p(:,:,i)=(eye(length(phi(i,:)))-k(:,i)*phi(i,:))*p(:,:,i-1);
    teta(:,i)=teta(:,i-1)+k(:,i)*(y(i)-phi(i,:)*teta(:,i-1));
    y_hat(i)=phi(i,:)*teta(:,i);
    E(i)=y(i)-y_hat(i);
end
%% validation 
Y_valid=y(150:200);
y_valid=y(150:200);
U_valid=U(150:200);
e_valid=e(150:200);
for j=3:50
      Y_valid(j)=-(A(2)*Y_valid(j-1)+A(3)*Y_valid(j-2))+B(1)*U_valid(j)+B(2)*U_valid(j-1)+C(1)*e_valid(j)+C(2)*e_valid(j-1)+C(3)*e_valid(j-2);
y_valid(j)=-(teta(1,i)*y_valid(j-1)+teta(2,i)*y_valid(j-2))+teta(3,i)*U_valid(j)+teta(4,i)*U_valid(j-1)+1*e(j)+teta(5,i)*e_valid(j-1)+teta(6,i)*e_valid(j-2);
ERROR_Valid(j)=abs(y_valid(j)-y(j));
end

subplot(2,3,1)
plot(3:i,A(2),'r*');grid on
hold on
plot(teta(1,:),'b--');legend('a1')
subplot(2,3,2)
plot(3:i,A(3),'r*');grid on
hold on
plot(teta(2,:),'b--');legend('a2')
subplot(2,3,3)
plot(3:i,B(1),'r*');grid on;grid on
hold on
plot(teta(3,:),'b--');legend('b1')
subplot(2,3,4)
plot(3:i,B(2),'r*');grid on
hold on
plot(teta(4,:),'b--');legend('b2')
subplot(2,3,5)
plot(3:i,C(1),'r*');grid on
hold on
plot(teta(5,:),'b--');legend('c1');grid on
subplot(2,3,6)
plot(3:i,C(2),'r*');grid on
hold on
plot(teta(6,:),'b--');legend('c2')
%

%% output train plot
figure;
plot(U,'linewidth',2);legend('INPUT');grid on
figure
plot(y,'linewidth',2);legend('OUTPUT');grid on
figure;
subplot(2,1,1)
plot(y,'linewidth',2);hold on
plot(y_hat,'r','linewidth',2);legend('sys ouput','y train');grid on
subplot(2,1,2)
plot(E,'r','linewidth',2);legend('Train ERROR');grid on

%% output validation 
figure;
subplot(2,1,1)
plot(Y_valid,'linewidth',2);hold on
plot(y_valid,'r','linewidth',2);legend('sys ouput','y validation');grid on
subplot(2,1,2)
plot(ERROR_Valid,'r','linewidth',2);legend('Validation ERROR ');grid on










