function [a]= gfun(alpha,x)
% x is a vector
%alpha is zero attracting parameter
[m,~]=size(x);
a=zeros(m,1);
for i=1:m
   if x(i)>1/alpha || x(i)<-1/alpha ||x(i)==0
       a(i)=0;       
   else 
       a(i)=1-alpha*abs(x(i));       
   end
end
            