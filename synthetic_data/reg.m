clc;clear;
N=512;% length of sparse presentation
L=300;mc=1;%  length of  signal
% SNR1=zeros(1,mc);SNR2=SNR1;
H=sqrt(1/L)*randn(L,N);% measurement matrix
% H=orth(H);
k =20;% nz component <N
sp=k/N; %sparsity degree
x = full(sprandn(N,1,sp)); %sparse presentation
x=x./norm(x);%normalize
alpha=1.5;gamma=5e-3;
n=starnd(alpha,gamma,0,0,L); %impulsive disturbance
s=H*x; %signal
y=s+n';% received signal
q1=0; q2=0; p=1.7;
% impulsive detection
imp_index=find(abs(y)>0.25);
y(imp_index)=0.25*sign(y(imp_index));
% ymax=max(abs(y));
% y=y./ymax;
 figure; plot(y,'g-');
%%
eta1=0:0.001:0.03;
eta2=0:0.002:0.06;
[lambda1,lambda2, MSD]=bestreg( eta1,eta2,y,q1,q2,p,H,x);
% [Lambda1,Lambda2]=meshgrid(lambda2,lambda1);
% figure
% mesh(Lambda1,Lambda2,MSD);
% figure;bar3(MSD);
lgMSD=10*log10(MSD);
figure;bar3(lgMSD);
figure; imagesc(lambda1, lambda2,lgMSD);
figure;contourf(lambda1,lambda2,lgMSD');