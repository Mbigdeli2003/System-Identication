function [U Y]=input_output
U=zeros(1,200);
n=200;
a=rand(1,n);
for i=1:200
    if a(i)<=.1
        U(i)=1;
    elseif a(i)<=.2 && a(i)>.1
        U(i)=5;
    elseif a(i)>.2 
        U(i)=3;
    end
    Y(i)=4*U(i).^2+7*U(i).^(1/2);
      
end

