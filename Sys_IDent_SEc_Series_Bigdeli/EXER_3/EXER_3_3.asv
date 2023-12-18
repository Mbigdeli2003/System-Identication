clc
clear all
close all
%% Training
%tarain number
N=150;
[U y e phi]=INPUT_OUTPUT;
INP=U(1:150)';
Y=y(1:150)';
U_Remover=eye(size(INP'*inv(INP*INP')*INP))-(INP'*inv(INP*INP')*INP);
G=(1/N)*Y*U_Remover*phi;
%% wwightening matrix to perform G
w2=inv(1/N*phi'*U_Remover*phi)*phi'*U_Remover;
w1=1;
G=w1*G*w2;
[u s v]=svd(G);
R=inv(phi*phi');
[n m]=find(s>0);
s1=s(n,m);
U1=u(:,1:max(m));
v1=v(:,1:max(m));
Or=U1*R;
C_hat=Or(1,:);
O_up=Or(1:max(size(Or))-1,:);
O_down=Or(2:max(size(Or)),:);
A_hat=inv(O_up'*O_down)*O_up'*O_down;
[Q R1]=qr([phi Y']);
R1=R1(1:min(size(phi)),1:min(size(phi)));
Y_hat=U1*R1*pinv(R1)*s1*v1;
X_hat=inv(Or)*Y_hat;
% B_hat=(X_hat-A_hat*X_hat)/INP;
syms q
INV_A=C_hat*inv(q*eye(size(A_hat))-A_hat);
[n,d]=numden(INV_A);
n_u=max(size(sym2poly(n(1,1))));
n_y=max(size(sym2poly(d(1,1))));
 phi_A=[-Y_hat(n_y-1) -Y_hat(n_y-2) -Y_hat(n_y-3) -Y_hat(n_y-4)  INP(150) INP(149) INP(148) INP(147)];
 YY=[Y_hat' 0 0 0 0];
KK=pinv(phi_A)*Y_hat(1,1);
B=KK(1:4);
D=KK(5:8);
YI=C_hat*X_hat+D(1,1)*INP;







% end


