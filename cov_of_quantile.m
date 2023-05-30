% 抽样-拟合-估计
% clear;clc;
close all
gev_para =   [-0.0408    ;1.2842   ;57.1802];
pd = makedist('gev','k',gev_para(1),'sigma',gev_para(2),'mu',gev_para(3));
% pd = makedist('gev','k',-0.1,'mu',3600,'sigma',1400);
p = logspace(-10,-1);
x = icdf(pd,1-p);
pdfx = pdf(pd,x);
y = p./(x.*pdfx);
figure;semilogx(p,y)
grid on
xlabel('$P_f$','Interpreter','latex')
ylabel('$\frac{P_f}{\xi f\left(\xi\right)}$','Interpreter','latex')
% pd = makedist('gev','k',0.1,'mu',3600,'sigma',1400);
% p = logspace(-10,-1);
% x = icdf(pd,1-p);
% pdfx = pdf(pd,x);
% y = p./(x.*pdfx);
% figure;semilogx(p,y)
% grid on
% xlabel('$P_f$','Interpreter','latex')
% ylabel('$\frac{P_f}{\xi f(\xi)}$','Interpreter','latex')
pd = makedist('gev','k',0,'mu',3600,'sigma',1400);
p = logspace(-10,-1);
x = icdf(pd,1-p);
pdfx = pdf(pd,x);
y = p./(x.*pdfx);
figure;semilogx(p,y)
grid on
xlabel('$P_f$','Interpreter','latex')
ylabel('$\frac{P_f}{\xi f(\xi)}$','Interpreter','latex')
pF = 4e-6;
p0 = 0.1;
m = ceil(log(pF)/log(p0));
ga = m;
N = 500;
% lower bound
r = 2;
delta = sqrt((1-p0)/p0*m^(r-1)*(1+ga)/N)
% upper bound
r = 3;
delta = sqrt((1-p0)/p0*m^(r-1)*(1+ga)/N)