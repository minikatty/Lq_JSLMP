function h=wiener_filter2(x,s,q)
% To calculate the coeficients of the Wiener filter 
% In case of missing values it ignores those samples
x=x(:);
s=s(:);
X=toeplitz(x,[x(1),zeros(1,q)]);
testgap=s'+sum(X');%to consider all NaNs in them
i=find(~isnan(testgap));
s=s(i);
X=X(i,:);
h=inv(X'*X)*(X'*s);

