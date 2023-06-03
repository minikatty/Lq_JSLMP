%%  alpha=1.5
%case1 gamma=1*1e-2;
ps=1:0.2:2;
case1=[4.845 4.802 5.185 5.867 5.852 5.568].*10^(-2);
case1=flip(case1);
%case2 gamma=4*1e-3;
case2=[8.92 9.15 10.14 10.64 10.58 10.30].*1e-3;
case2=flip(case2);
%case3 gamma=4*1e-3;
case3=[1.090 0.8466 0.5644  0.7856  0.7472 1.265].*1e-3;
case3=flip(case3);
figure;
semilogy(ps,case1);
hold on;
semilogy(ps,case2);
hold on;
semilogy(ps,case3);

%%
clc;clear;
ps=1:0.2:2;
case1=[8.61 8.466 8.644 8.782 8.531 8.668].*10^(-2);
case1=flip(case1);
%case2 gamma=4*1e-3;
case2=[1.72 1.738 1.682 1.676 1.707 1.719].*1e-2;
case2=flip(case2);
%case3 gamma=4*1e-3;
case3=[7.269 6.457 3.3449 3.0500  3.252 1.999].*1e-4;
case3=flip(case3);
figure;
semilogy(ps,case1);
hold on;
semilogy(ps,case2);
hold on;
semilogy(ps,case3);