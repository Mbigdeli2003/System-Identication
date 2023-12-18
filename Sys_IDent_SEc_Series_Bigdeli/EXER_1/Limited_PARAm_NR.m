clc
clear all
close all
%% initial condition no#4
% function calling
[U Y]=input_output;

%%%%%%%%%%% training %%%%%%%%%%%%%%%%
n=150;
init=[2 1.7 5 1]; % initial condition
teta1=init(1);
teta2=init(2);
teta3=init(3);
teta4=init(4);
teta=[teta1;teta2;teta3;teta4];
alpha=.4*ones(1,n);

for i=1:n
    
 df=[U(i)^teta(2,i);teta(1,i)*log(U(i))*U(i)^teta(2,i);U(i)^teta(4,i);teta(3,i)*log(U(i))*U(i)^teta(4,i)]*(teta(1,i)*U(i)^teta(2,i)+teta(3,i)*U(i)^teta(4,i)-Y(i));
    d2f=[0 log(U(i))*U(i)^teta(2,i) 0 0;log(U(i))*U(i)^teta(2,i) teta(1,i)*log(U(i)^2*U(i)^teta(2,i)) 0 0;0 0 0 log(U(i))*U(i)^teta(4,i);0 0 log(U(i))*U(i)^teta(4,i) teta(3,i)*(log(U(i))^2*U(i)^(teta(4,i)))]*(teta(1,i)*U(i)^teta(2,i)+teta(3,i)*U(i)^teta(4,i)-Y(i))+...
        [U(i)^teta(2,i);teta(1,i)*log(U(i))*U(i)^teta(2,i);U(i)^teta(4,i);teta(3,i)*log(U(i))*U(i)^teta(3,i)]*[U(i)^teta(2,i) teta(1,i)*log(U(i))*U(i)^(teta(2,i)) U(i)^teta(4,i) teta(3,i)*log(U(i))*U(i)^(teta(4,i))];
    teta(:,i+1)=teta(:,i)-(alpha(i)*pinv(d2f')*df);
    if teta(2,i+1)<1
        teta(2,i+1)=1;
    elseif teta(2,i+1)>2
        teta(2,i+1)=2;
    end
       if teta(4,i+1)<0
        teta(4,i+1)=0;
    elseif teta(4,i+1)>1
        teta(4,i+1)=1;
    end
        

Y_train(i)=teta(1,i)*U(i)^teta(2,i)+teta(3,i)*U(i)^teta(4,i);
   ERR_train(i)=abs(Y_train(i)-Y(i));
end
SSE_train=sum((Y_train-Y(1:150)).^2)

%% validation
U_VALID=U(151:200);
Y_sys=Y(151:200);
for j=1:50
Y_valid(j)=teta(1,i)*U_VALID(j)^teta(2,i)+teta(3,i)*U_VALID(j)^teta(4,i);
ERR(j)=abs(Y_valid(j)-Y_sys(j));

end
SSE=sum((Y_valid-Y_sys).^2)
%% plotting
% Input plot
figure;
plot(U,'linewidth',2.5);grid on
%% parameter Converegence

%% Training Plot
 figure;
 subplot(2,1,1);
 plot(ERR_train,'b','linewidth',2.5);legend('train error ');
 subplot(2,1,2);
 plot(Y_train,'b--','linewidth',3);hold on
 plot(Y(1:150),'r--','linewidth',2.25);legend('Y train','Y')

    

%% validation Plot
 figure;
 subplot(2,1,1);
 plot(ERR,'b','linewidth',2.5);legend('validation error ');grid on
 subplot(2,1,2);
 plot(Y_valid,'b--','linewidth',3);hold on
 plot(Y_sys,'r--','linewidth',2.25);legend('Y Valid','Y');grid on











