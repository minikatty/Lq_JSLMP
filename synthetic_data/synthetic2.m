%% signal model: y=Hx+n
clc;clear;
N=512;% length of sparse presentation
L=300;mc=10;%  length of  signal
SNR1=zeros(1,mc);SNR2=SNR1;
min=0;
lambda1=0:0.0025:0.025;
lambda2=0:0.005:0.05;
i=0;j=0;
p=1.7;
%% Monte Carlo
for k=5:5:45   
            i=i+1;              
            for l=1:mc                    
                    H=sqrt(1/L)*randn(L,N);% measurement matrix                                                            
                    sp=k/N; %sparsity degree
                    x = full(sprandn(N,1,sp)); %sparse presentation
                    x=x./norm(x);%normalize
                    alpha=1.5;gamma=4e-3;
                    n=starnd(alpha,gamma,0,0,L); %impulsive disturbance
                    s=H*x; %signal
                    y=s+n';% received signal
%%  Adam+VSS JSLMP
                    for q1=0:0.2:1
                        j=j+1;
                        q2=q1;
                        [~,~,~,MSD(j,l)]=bestreg(y,q1,q2,p,H,x,lambda1,lambda2);  %best performance                                                               
                    end
                    j=0;
            end
             count=sum(MSD<0.01,2);
             rate(:,i)=count./mc;
end
            figure;plot(5:5:45,rate(1,:));
