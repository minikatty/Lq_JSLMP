clc;clear;
mc=2000; %(mc>1)
Mean=1/mc;mu=0.0000001;lambda=[0.006 0.02 0.04 0.04  0.06 0.075]';
cNZ=30;M=300;N=512;
LL=M+M-cNZ;
ps=1:0.2:2;np=length(ps);sigma=1;
M1=zeros(np,M);M2=M1;M3=M2;
alpha=1.5;gamma=4e-3;
for i=1:np
p1=2*ps(i)-2;
p2=2*ps(i)-4;
p3=ps(i)-2;
for j=1:mc
n=starnd(alpha,gamma,0,0,M);
M1(i,:)=Mean*abs(n).^p1+M1(i);
M2(i,:)=Mean*abs(n).^p2+M2(i);
M3(i,:)=Mean*abs(n).^p3+M3(i);
end
end

means1=mean(M1,2);
means2=mean(M2,2);
means3=mean(M3,2);
C3=1-mu*sigma^2.*means2;
C4=2-(2+cNZ)*mu*(sigma)^2.*means2;
si=(C4.*LL*mu^2*sigma^2.*means1+cNZ*mu^3*sigma^4.*means1.*means2)./...
    (2*mu*sigma^2.*C3.*C4-2.*C3*LL*mu^2*sigma^4.*means2+sqrt(8/pi)*mu.*lambda.*...
    C3.*C4*LL);
kappa=mu.*lambda;
% term4=2*C2*kappa.^2./(C4*mu^2*sigma^4.*means3);
theory=2.*C3*LL.*si./C4+cNZ*mu.*means1./C4+cNZ*kappa.^2./(C4*mu*sigma^2);
figure;plot(ps,theory);
% figure;semilogy(ps,means1);figure;semilogy(ps,means2);
% figure;semilogy(ps,means3);