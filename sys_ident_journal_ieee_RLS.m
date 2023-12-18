%% Hierarchical Least Squares Identification for Linear SISO Systems With Dual-Rate Sampled-Data
%% morteza bigdeli 9023405064_journal paper implement 
%% ''RLS''
clc
clear all
close all
format short
syms z
%% input defining
n=12000;
var_n=(1^(1/2));%  variance
U=random('normal',0,var_n,n,1);
varN=(0.25)^(1/2);
v=random('normal',0,varN,n,1);

%% Defining system parameters
% A(z^-1)=1-1.1z^-1+0.5z^-2-0.37z^-3
a1=-1.1;a2=0.5;a3=-0.37; 
A=[1 a1 a2 a3];
AA=flipdim(A,2);
%B(z^-1)=0.46z^-1-0.34z^-2+0.5z^-3
b0=0;b1=0.46;b2=-0.34;b3=0.5;
B=[b0 b1 b2 b3];
BB=flipdim(B,2);
%D(z^-1)=1-a1(z^-1)+a2(z^-2)-a3(z^-3);
d1=a1;d2=-a2;d3=a3;
%% define delta_ns (var(W)/var(X))
w=idpoly(A,1);
W=sim(w,v);
x=idpoly(A,B);
X=sim(x,U);
delta_ns=((var(W)/var(X))^(1/2))*100



%% polynominals Trans 
%% Creating alfa & Beta taking q=2
% A(z^-1)=1-1.1z^-1+0.5z^-2-0.37z^-3
AZ=poly2sym(AA,z^-1);
%B(z^-1)=0.46z^-1-0.34z^-2+0.5z^-3
BZ=poly2sym(BB,z^-1);
%D(z^-1)=(1+z1z^-1)(1+z2z^-1)(1+z3z^-1);z1,z2,z3 roots of A(z^-1);
Z=roots(A);
DZ=(1+Z(1,1)*z^(-1))*(1+Z(2,1)*(z^-1))*(1+Z(3,1)*(z^-1));
DZ=expand(DZ);
DZ=vpa(DZ);
% alfa(z^-1)=D(z^-1)A(z^-1)
alfaZ=DZ*AZ;
alfaZ=expand(alfaZ);
alfaZ=vpa(alfaZ);
% % % disp('alfa(z^-1)')
% % % pretty(alfaZ)
% beta(z^-1)=D(z^-1)B(z^-1)
betaZ=DZ*BZ;
betaZ=expand(betaZ);
betaZ=vpa(betaZ);
% % % disp('beta(z^-1)')
% % % pretty(betaZ)
%% determining coeffs
% alfa coeffs
Alfa_Z=z^6*alfaZ;
Alfa=sym2poly(Alfa_Z);
% beta coeffs
Beta_Z=z^6*(betaZ);
Beta=sym2poly(Beta_Z);
% D coeffs
DZ=z^3*DZ;
D=sym2poly(DZ);
%% RLS with kq
q=2;
na=length(A)-1;
nb=length(B);
nd=(length(D)-1);
phi=zeros(14,13);
alpha =1e6;
P(:,:,12) = alpha*eye(length(phi(1,:)));
yk=zeros(1,length(U));
nv_k_hat=repmat((1/alpha)*ones(13,1),1,14);
sa_hat=zeros(14,10);
teta=[A(2:end) B]';
teta_hat=(1/alpha)*ones(7,13);
delta=zeros(1,length(U));
%% RLS
for k=7:(length(U)/2)-2
     for l=1:na
    
    yk(k*q)=-(Alfa(1,3)*yk((k*q)-2)+Alfa(1,5)*yk((k*q)-4)+Alfa(1,7)*yk((k*q)-6))+Beta(1,1)*U((k*q)-1)+Beta(1,2)*U((k*q)-2)+Beta(1,3)*U((k*q)-3)+Beta(1,4)*U((k*q)-4)+Beta(1,5)*U((k*q)-5)+Beta(1,6)*U((k*q)-6)+D(1,1)*v(k*q)+D(1,2)*v((k*q)-1)+D(1,3)*v((k*q)-2)+D(1,4)*v((k*q)-3);
  
  yk(k*q+1) = interp1([k*q k*q+q],[yk(k*q) yk(k*q+q)],k*q+1);

  phi(k*q,:)=[-yk((k*q)-q) -yk(k*q-(2*q)) -yk(k*q-(3*q)) U(k*q) U(k*q-1) U(k*q-2) U(k*q-3) U(k*q-4) U(k*q-5) U(k*q-6) v(k*q-1) v(k*q-2) v(k*q-3)];    
    P(:,:,k*q)=P(:,:,k*q-q)-((P(:,:,k*q-q)*phi(k*q,:)'*phi(k*q,:)*P(:,:,k*q-q))/(1+phi(k*q,:)*P(:,:,k*q-q)*phi(k*q,:)'));



  sa_hat(k*q-l,:)=[-yk(k*q-l-1) -yk(k*q-l-2*q) -yk(k*q-l-3*q) U(k*q-l) U(k*q-l-1) U(k*q-l-2) U(k*q-l-3) U(k*q-l-4) U(k*q-l-5) U(k*q-l-6) ]; 
 phi(k*q-l,:)=[sa_hat(k*q-l,:) v(k*q-l-1) v(k*q-l-2) v(k*q-l-nd)];
   v(k*q-l)=yk(k*q-l)-phi(k*q-l,:)*nv_k_hat(:,k*q-q);


 nv_k_hat(:,k*q)=nv_k_hat(:,k*q-q)+P(:,:,k*q)*phi(k*q,:)'*(yk(k*q)-phi(k*q,:)*nv_k_hat(:,k*q-q));
  nv_k_hat(:,(k*q)+1)=nv_k_hat(:,k*q);
S_a(:,:,k*q)=[1                                    0                  0                 0      ;
              0                                    1                  0                 0     ;
              nv_k_hat(1,k*q)                      0                  1                 0     ;
              0                                    nv_k_hat(1,k*q)    0                 1      ;
              nv_k_hat(2,k*q)                      0                  nv_k_hat(1,k*q)   0       ;
              0                                    nv_k_hat(2,k*q)    0                 nv_k_hat(1,k*q);
              nv_k_hat(3,k*q)                      0                  nv_k_hat(2,k*q)     0;
              0                                    nv_k_hat(3,k*q)    0                 nv_k_hat(2,k*q);
              0                                    0                  nv_k_hat(3,k*q)   0;
              0                                    0                  0                 nv_k_hat(3,k*q)];

  
  
  
  S_b(:,:,k*q)=[0                      0                   0;
    nv_k_hat(4,k*q)                    0                   0;
    nv_k_hat(5,k*q)                    nv_k_hat(4,k*q)     0;
    nv_k_hat(6,k*q)                    nv_k_hat(5,k*q)     nv_k_hat(4,k*q);
    nv_k_hat(7,k*q)                    nv_k_hat(6,k*q)     nv_k_hat(5,k*q);
    nv_k_hat(8,k*q)                    nv_k_hat(7,k*q)     nv_k_hat(6,k*q);
    nv_k_hat(9,k*q)                    nv_k_hat(8,k*q)     nv_k_hat(7,k*q);
    nv_k_hat(10,k*q)                   nv_k_hat(9,k*q)     nv_k_hat(8,k*q);
    0                                  nv_k_hat(10,k*q)    nv_k_hat(9,k*q);
    0                                  0                   nv_k_hat(10,k*q)];
pro(:,k*q)=[nv_k_hat(4:10,k*q)' 0 0 0];
S(:,:,k*q)=[-S_b(:,:,k*q) S_a(:,:,k*q)];


teta_hat(:,k*q)=pinv(S(:,:,k*q))*pro(:,k*q);
teta_hat(:,k*q+1)=teta_hat(:,k*q);
delta(k)=norm((teta_hat(:,k)-teta),1)/norm(teta,1);
end 
     
end
%% Plot
figure;
plot(delta,'linewidth',2);legend('delta');axis([0 k 0 0.7]);grid on;

figure;
subplot(2,3,1)
plot(1:2*k+1,A(2),'r*');
hold on
plot(teta_hat(1,:),'b--');legend('a1')
subplot(2,3,2)
plot(1:2*k+1,A(3),'r*')
hold on
plot(teta_hat(2,:),'b--');legend('a2')
subplot(2,3,3)
plot(1:2*k+1,A(4),'r*')
hold on
plot(teta_hat(3,:),'b--');legend('a3')
subplot(2,3,4)
plot(1:2*k+1,B(2),'r*')
hold on
plot(teta_hat(5,:),'b--');legend('b1')
subplot(2,3,5)
plot(1:2*k+1,B(3),'r*')
hold on
plot(teta_hat(6,:),'b--');legend('b2')
subplot(2,3,6)
plot(1:2*k+1,B(4),'r*')
hold on
plot(teta_hat(7,:),'b--');legend('b3')




















