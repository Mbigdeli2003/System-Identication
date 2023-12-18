function [u y]=input_output
close all; clear; clc;
u=zeros(1,200);
n=200;
a=rand(1,n);
for i=1:200
    if a(i)<=.1
        u(i)=1;
    elseif a(i)<=.2 && a(i)>.1
        u(i)=5;
    elseif a(i)>.2 
        u(i)=3;
    end
    y(i)=4*u(i).^2+7*u(i).^(1/2);
      
end
