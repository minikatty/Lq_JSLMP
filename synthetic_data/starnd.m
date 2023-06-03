function Y=starnd(alpha,gamma,beta,a,N) 


% STARND draws from standard stable distributions 
% (Chambers-Mallow-Stuck method). 
% 
% X=starnd(alpha,beta,N) 
% 
% Inputs: 
% - alpha is the characteristic exponent (0<alpha<=2) 
% - beta is skewness parameter (-1<beta<1) 
% - N is the number of samples 
% 
% Output: 
% - X is a 1 x N vector containing the sampled values 
% 
% Linked m-files: this function uses the gamrnd.m Statistics Matlab Toolbox 
% function for sampling from Gamma distribution. 
% 
% Author: C. FEVOTTE cf269 AT cam DOT ac DOT uk 
% (please report any bug) 
% Reference: R. Weron, "On the Chambers-Mallows-Stuck method for 
% simulating skewed stable random variables", Technical Report, Hugo 
% Steinhaus Center for Stochastic Methods, Wroclaw, 1996 
if nargin==4; 
    N=1; 
end 

V=rand(1,N)*pi-pi/2; 
W=gamrnd(ones(1,N),1); 


if alpha==1 
    X=(2/pi)*((pi/2+beta*V).*tan(V)-beta*log((pi/2)*W.*cos(V)./(pi/2+beta*V)));
    Y=gamma*X+(2/pi)*beta*gamma*log(gamma)+a;
else 
   B=atan(beta*tan(pi*alpha/2))/alpha; 
   S=(1+beta^2*tan(pi*alpha/2)^2)^(1/(2*alpha)); 
   X=S.*sin(alpha*(V+B))./(cos(V).^(1/alpha)).*(cos(V-alpha*(V+B))./W).^((1-alpha)/alpha); 
   Y=gamma*X+a;
end 