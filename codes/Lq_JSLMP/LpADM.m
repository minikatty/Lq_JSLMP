function [res] =LpADM(y,lamda,p,x0,rho)
% l1_lp_admm solves (without smoothing)
%
%   minimize || Ax - y ||_p^p + \lamda || x ||_1
%
% Inputs:
%	A: DCT matrix
%	y: speech frame
%	lamda: regularization parameter 
%	x0: initialization 
% Outputs
%	x: the recovery
%Convergence setup
L=length(y);m=L;n=L;
A=dctmtx(L);
Max=max(abs(y));
if Max==0
    Max=1e-8;
end
y=y./Max;
max_iter = 200;
ABSTOL = 1e-6;
if nargin<5
	rho = 1;
end    
%Initialize
if nargin<4
	x  = zeros(n,1);
else
	x = x0;
end

v = zeros(m,1);
w = zeros(m,1); 
% xm1 = x;
rhoT = rho;
rho = 10;

% out.e  = [];
% out.et = [];out.f = [];

for i = 1 : max_iter
    vm1 = v;    
    if rho<rhoT
        rho =rho*1.03;
    end           
    %v-step
	tv = A*x-y-w/rho;
    v  = shrinkage_Lp(tv, p, 1/lamda, rho);    
    %x-step
    tao = 0.9; % for orthonornal A
    z = x - tao*(A'*(A*x-y-v-w/rho)); 
    x = sign(z) .* max(abs(z)-tao/rho,0);    
	%w-step
    Ax = A*x;
	w = w - rho*(Ax - y - v);    
%     xm1 = x;
    %terminate when both primal and dual residuals are small
    if (norm(rho*(v-vm1))< sqrt(n)*ABSTOL && norm(Ax-y-v)< sqrt(n)*ABSTOL) 
        break;
    end
end
res=A*x.*Max;
end
