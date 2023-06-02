%*******************************************************************************
%                                                                                                                             *
%                                          Conbined ADAM+VSS                                             *
%                                                                                                                             *
%********************************************************************************
function [VSSLMS_ADM,index,nv]=VSSLMS_Adam_l0(H,y,iter_max,tol,x_true)% l0  approximate 
[L,N]=size(H);
%initial condition
n=1;mu=0.4;mu_max=0.4;mu_min=1e-5;
lambda=2e-5; % condition to end ADM
tol2=3.5e-5;alpha2=3;
w_old=zeros(N,1);%w_new=zeros(N,1);
m_old=zeros(N,1);w_new=m_old;
v_old=zeros(N,1);
alpha=1e-3;beta1=0.992;beta2=1-1e-3;epsilon=1e-8;

%---------------------------------Adam------------------------------------
%****************************************************************************
    while (n~=1 && norm(w_new-w_old)>tol && n<=iter_max) || n==1
    k=mod(n,L)+1; 
    d=y(k);                                          % received signal
    X=H(k,:)';
    w_old=w_new;
     
    g=-lambda*g_function(10,w_old)-X*(d-X'*w_old);   % gradient    
    m_new=beta1*m_old+(1-beta1)*g;                            % update bias first moment estimate    
    v_new=beta2*v_old+(1-beta2)*(g.*g);                        % update biased second raw moment estimate
    m1=m_new/(1-beta1^(n+1));                                      % compute biased-corrected first moment estimate 
    v1=v_new/(1-beta2^(n+1));                                         % compute biased-corrected second raw moment estimate
    nv(n)=norm(v1);
%    alpha1=alpha*sqrt(1-beta2^(n+1))/(1-beta1^(n+1));% new stepsize                   
    w_new=w_old-alpha*m1./(sqrt(v1)+epsilon);            % update parameters 
    % for debug -alpha1*m1./(sqrt(v1)+epsilon)
     VSSLMS_ADM.MSD(n)=sum((x_true-w_new).^2);
    m_old=m_new;
    v_old=v_new;       
    if nv(n)<=tol2&&(n>=10)
        index=n; 
       break;        
    end 
    n=n+1;  
    end
    n=n+1; 
%---------------------------------------ADAM-END-------------------------------------

%----------------------------------------VSSLMS----------------------------------------
while ( n<=iter_max) %norm(w_new-w_old)>tol &&
    w_old=w_new;
    mu_old=mu;
    k=mod(n,L)+1;      %avoid deviding by 0
    d=y(k);
    X=H(k,:)';
    e=d-X'*w_old;
    %-------------------------------------------------
    % VSS mechanism   
    mu_new=0.999*mu_old+4e-3*e^2;
%     mu_new=0.95*mu_old+5e1*e^2;
    if  mu_new>mu_max
        mu=mu_max;
    elseif mu_new<mu_min
        mu=mu_min;
    else
        mu=mu_new;    
    end    
    gamma=mu*lambda;
    %-------------------------------------------------
    w_new=w_old+mu*e*X;         %gama= mu*lambda
    w_new=w_new+gamma*g_function(alpha2,w_old);
   VSSLMS_ADM.MSD(n)=(norm(x_true-w_new))^2;
    n=n+1;    
end
VSSLMS_ADM.x=w_new;


  
