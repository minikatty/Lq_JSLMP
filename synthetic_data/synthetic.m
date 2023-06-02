%% signal model: y=Hx+n
clc;clear;
N=512;% length of sparse presentation
L=300;mc=1;%  length of  signal
SNR1=zeros(1,mc);SNR2=SNR1;
%% Monte Carlo
for i=1:mc
H=sqrt(1/L)*randn(L,N);% measurement matrix
% H=orth(H);
k =40;% nz component <N
sp=k/N; %sparsity degree
x = full(sprandn(N,1,sp)); %sparse presentation
x=x./norm(x);%normalize
alpha=1.5;gamma=4e-3;
n=starnd(alpha,gamma,0,0,L); %impulsive disturbance
s=H*x; %signal
y=s+n';% received signal
% SNR1(i)=10*log10((sum(s.^2))/(2*(gamma^(2/alpha)))); % SNR defined in \alpha stable theory
% SNR2(i)=10*log10((sum(s.^2))/sum(n.^2));  % traditional SNR  definition
% SNR3=10*log10(sum(s.^2)/sum(n.^2));
%  figure; plot(y,'b-');
%  figure;plot(s,'r.-');
%% impulsive detection
imp_index=find(abs(y)>0.25);
y(imp_index)=0.25*sign(y(imp_index));
% ymax=max(abs(y));
% y=y./ymax;
%  figure; plot(y,'g-');
%%  Adam+VSS JSLMP
q1=0; q2=q1; p=1.5;
lambda1=0.03;lambda2=lambda1;
% res1=JSLMSP(y,q1,q2,p,H,x,lambda1,lambda2,'normal');
res2=JSLMSP3(y,q1,q2,p,H,x,lambda1,lambda2,'normal');
% initial.x=res1.para;
% initial.v=res1.v;
% res2=JSLMSP3(y,q1,q2,p,H,x,lambda1,lambda2,'normal',initial);
% figure;plot(res2.msd);
 if i==1
%         MSD1=1/mc.*res1.msd;
        MSD2=1/mc.*res2.msd;
 else
%         MSD=1/mc.*res1.msd+MSD;
        MSD2=1/mc.*res2.msd+MSD;
 end
% y_re=H*res.para*ymax;
%  figure;
%  semilogy(res1.msd,'r.-');hold on;
%  semilogy(10^(-2)*ones(length(MSD2),1),'g.-');
%  hold on;semilogy(res2.msd,'b-');
%   figure;semilogy(MSD,'r-');
end
figure;semilogy(MSD2,'r-');
% hold on;semilogy(MSD2,'b-');
hold on;semilogy(10^(-2)*ones(length(MSD2),1),'g.-');
disp(['change happened at:  ', num2str(res2.label)]);
disp(['VSS is used or not: ', num2str(res2.flag)]);
% fprintf('epsilonr is : %f\n', res.epsilonr);
 
