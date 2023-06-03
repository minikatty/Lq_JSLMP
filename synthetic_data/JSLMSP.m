%--------------------------------------------------------
%objective function:
%     J(w(n),v(n))=1/p*[y(n)-x(n)^T*w(n)-v(n)_k]^p+
%   lambda_1||w(n)||_q1^q1+lambda_2||v(n)||_q2^q2
%                 0<=q2,q2<=1
%                      1<p<2  
%--------------------------------------------------------
function [JSLMP]=JSLMSP(y,q1,q2,p,H,w_true,lambda1,lambda2,sym,initial)
   
L=size(H,1);iter_max=8e4;tol=1e-8;N=size(H,2);
n=1;flag=1;I=eye(L);
MSD=zeros(iter_max,1);thres=MSD;
lambda1=lambda1/L;lambda2=lambda2/L;        %adaptive filter parameters
flag2=strcmp(sym,'reweighted');

if  nargin<10
    w_prev=zeros(N,1); w_new=w_prev; v_prev=zeros(L,1);               %primary adaptive filter tap weights
else                                                 
    w_prev=initial.x;w_new=w_prev; v_prev=initial.v;
end 
                                                %secondary adaptive filter tap weights
%H=dctmtx(L);                                                       %speech deniose  with DCT
    tau=1e-8;                                                         % cliping threshold 
%------------------------Adam Parameter+ VSS Setting---------------------------
m1=zeros(N,1);v1=m1;m2=zeros(L,1);v2=m2;epsilon1=1e-8;
alpha=1e-3;beta1=0.9;beta2=1-1e-2;epsilon=1e-8;alpha1=10;
mu=0.01;mu_max=mu;mu_min=1e-5;mu_pre=mu; label=0;epsilonr=1;
%-----------------------------------------------------------------------------
if (q1==0)
    g1=@(X,e,w_prev)(-X*abs(e)^(p-1)*sign(e)+lambda1*alpha1*...
    sign(w_prev).*gfun(alpha1,w_prev));
elseif (q1==1)
    g1=@(X,e,w_prev)(-X*abs(e)^(p-1)*sign(e)+lambda1.*sign(w_prev));    
else
    if  flag2
        D=[H I];
        init=pinv(D)*y;  %least l2-norm solution   
        w_prev=init(1:N);  
        w_new=w_prev; %weiths initial
        w1=@(w1,q1,epsilonr)(w1.^2+epsilonr).^(q1/2-1);       % weights  matrix update function
        weight1=ones(N,1);
        g1=@(X,e,w_prev,W1)(2*lambda1*W1.*w_prev-abs(e)^(p-1)*sign(e)*X);
    elseif strcmp(sym,'normal')
    g1=@(X,e,w_prev)((lambda1*q1./(abs(w_prev).^(1-q1)+epsilon1))...
    .*sign(w_prev)-X*abs(e)^(p-1)*sign(e)); 
    else        
       error('argument ''sym''  is wrong.');
    end
end

if (q2==0)
     g2=@(Ik,e,v_prev)(-Ik*abs(e)^(p-1)*sign(e)+lambda2*alpha1*...
    sign(v_prev).*gfun(alpha1,v_prev));
elseif (q2==1)
    g2=@(Ik,e,v_prev)(-Ik*abs(e)^(p-1)*sign(e)+lambda2.*sign(v_prev));
else
    if flag2
        v_prev=init(N+1:end);                %weiths initial
        w2=@(w2,q2,epsilonr)(w2.^2+epsilonr).^(q2/2-1);       % weights  matrix update function
        weight2=ones(L,1);
        g2=@(Ik,e,v_prev,W2)(2*lambda2*W2.*v_prev-Ik*abs(e)^(p-1)*sign(e));
    elseif strcmp(sym,'normal')
        g2=@(Ik,e,v_prev)((lambda2*q2./(abs(v_prev).^(1-q2)+epsilon1))...
    .*sign(v_prev)-Ik*abs(e)^(p-1)*sign(e));
    else
        error('argument ''sym''  is wrong.');
    end
end
v_hat_prev1=zeros(N,1);v_hat_prev2=zeros(L,1);
while (n~=1 && norm(w_prev-w_new)>tol && n<=iter_max) || n==1     
    k=mod(n,L)+1; 
    d=y(k);
    X=H(k,:)'; Ik=I(k,:)';    
    w_prev=w_new;
    MSD(n)=sum((w_prev-w_true).^2);
    e=d-X'*w_prev-Ik'*v_prev;
    if ~flag2
    gra1=g1(X,e,w_prev); % primary gradient
    gra2=g2(Ik,e,v_prev); % secondary gradient  
    else
    gra1=g1(X,e,w_prev,weight1); % primary gradient
    gra2=g2(Ik,e,v_prev,weight2); % secondary gradient  
    end
    if flag                
        m1=beta1*m1+(1-beta1)*gra1;                  % update bias first moment estimate
        v1=beta2*v1+(1-beta2)*(gra1.*gra1);            % update biased second raw moment estimate
        m11=m1/(1-beta1^(n+1));                  % compute biased-corrected first moment estimate
        v11=v1/(1-beta2^(n+1));                     % compute biased-corrected second raw moment estimate
        %alpha1_new(n)=alpha*sqrt(1-beta2^(n+1))/(1-beta1^(n+1));
        v_hat1 = max(v_hat_prev1, v11);
        w_new=w_prev-alpha*m11./(sqrt(v_hat1)+epsilon);% update parameters 
%         w_new=w_prev-alpha*m11./(sqrt(v11)+epsilon);% update parameters 
        v_hat_prev1=v_hat1;                              
        m2=beta1*m2+(1-beta1)*gra2;                  % update bias first moment estimate
        v2=beta2*v2+(1-beta2)*(gra2.*gra2);         % update biased second raw moment estimate
        m21=m2/(1-beta1^(n+1));                          % compute biased-corrected first moment estimate
        v21=v2/(1-beta2^(n+1));
        v_hat2 = max(v_hat_prev2, v21);
        v_new=v_prev-alpha*m21./(sqrt(v_hat2)+epsilon);     
%         v_new=v_prev-alpha*m21./(sqrt(v21)+epsilon);    
        v_prev=v_new; 
        v_hat_prev2=v_hat2;
        thres(n)=sum(m11.^2);        
         if (thres(n)<=tau) && (n>3e3)
                flag=0;
                label=n;
                continue;
         end         
    else        
%-------------------------------------------------
    % VSS mechanism      
        mu_new=0.999*mu_pre+2e-2*e^2;
%     mu_new=0.95*mu_old+5e1*e^2;
         if  mu_new>mu_max
                mu=mu_max;
        elseif mu_new<mu_min
                mu=mu_min;
         else
                mu=mu_new;    
         end
         mu_pre=mu_new;
    %-------------------------------------------------     
            w_new=w_prev-mu*gra1;       
            v_new=v_prev-mu*gra2;
            v_prev=v_new;
    end
    if flag2
        if  norm(w_new-w_prev)<=sqrt(epsilonr)/100 && epsilonr>10^(-1)
            epsilonr=epsilonr/10;
        end
        weight1=w1(weight1,q1,epsilonr);
        weight2=w2(weight2,q2,epsilonr);
    end 
    n=n+1;    
end    
JSLMP.para=w_new;
JSLMP.v=v_new;
JSLMP.msd=MSD;
JSLMP.flag=flag;
JSLMP.label=label;
JSLMP.thres=thres;
JSLMP.epsilonr=epsilonr;
end
