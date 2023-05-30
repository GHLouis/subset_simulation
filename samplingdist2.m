% 抽样-拟合-估计
clear;clc;close all
pd = makedist('normal','mu',0,'sigma',1);
name = pd.DistributionName;
pf_ex = 1 - exp(log(1-1/75)/(250));
q = 1 - pf_ex;
% confidence interval
alpha = 1 - 0.95;
k = abs(norminv(alpha/2));
% 误差上限
EL = 0.01;
fprintf('具有%g%%保证率下相对误差不超过%g%%的分位值估计 \n', 100*(1-alpha), 100*EL)
%% 理论值
Q_mean = icdf(pd,q);
% N = (k/EL)^2 * (sqrt(q*(1-q))/(pdf(pd, Q_mean) * Q_mean))^2;
% N = ceil(N);
% Q_std = sqrt(q*(1-q)/N)/pdf(pd, Q_mean);
% Q_cov = Q_std/Q_mean;
% fprintf('# 理论计算 \n')
% fprintf('分位值为 %g \n', Q_mean)
% fprintf('理论上变异系数为 %g \n', Q_cov)
% fprintf('需要的总样本数为 %g \n', N)
%% 基于Bootstrap方法估计分位值的标准差
% Bootstrap方法估计分位值
n = 2000;
x = random(pd,[n, 1]);
B = 1000; % bootstrap样本数量
q_hat = zeros(B, 1); % 存储所有分位数的向量
for i = 1:B
  boot_sample   = datasample(x, n, 'Replace', true);
  pd_fit        = fitdist(boot_sample, name);
  q_hat(i)      = icdf(pd_fit, q);
end
% 估计值的变异系数
c = std(q_hat)/mean(q_hat);
% 根据均值抽样分布，确定需要重复的次数
N = max(ceil((k/EL*c)^2) ,1);
fprintf('# Bootstrap \n')
fprintf('给定%g个样本 \n', n)
% fprintf('理论上变异系数为%g \n',Q_cov)
fprintf('Bootstrap方法估计的变异系数为 %g \n', c)
fprintf('因此需要重复%g次 \n', N)
fprintf('需要的总样本数为 %g \n', n*N)
%% result
x       = random(pd,[n*N, 1]);
q_hat = zeros(1,1);
for i = 1:1
   pd_fit       = fitdist(x(:,i), name);
  q_hat(i)      = icdf(pd_fit, q); 
end
fprintf('# Single run \n')
Q_es_mean   = mean(q_hat);
Q_es_std    = std(q_hat);
Q_es_cov    = std(q_hat)/mean(q_hat);
er_100      = (Q_es_mean - Q_mean)/Q_mean*100;
fprintf('使用的样本数为 %g \n', n*N)
fprintf('Estimated quantile is %g\n', Q_es_mean)
fprintf('Relative bias %g%% \n', er_100)

%% 95% confidence interval
Np = 20;
x = random(pd,[1000, Np]);
% fprintf('# *** Multiple runs ***\n')
q_hat = zeros(Np,1);
for i = 1:Np
   pd_fit       = fitdist(x(:,i), name);
  q_hat(i)      = icdf(pd_fit, q);
end
error_ci_plot(q_hat,Q_mean)
% e       = (q_hat - Q_mean)/Q_mean * 100;
% dist    = histfit_pdf(e(:),name);
% xlabel('Relative bias (%)')
% ylabel('PDF')
% e_mean  = mean(e);
% e_std   = std(e);
% e_cov   = e_mean/e_std;
% % er_100  = (Q_es_mean - Q_mean)/Q_mean*100;
% fprintf('Mean of quantile is %g \n', mean(q_hat))
% fprintf('Mean of relative bias is %g%% \n', e_mean)
% fprintf('COV of relative bias is %g \n', e_cov)
% 
% ax      = gca;
% e_95L   = mean(e) - 2*std(e);
% e_95U   = mean(e) + 2*std(e);
% fprintf('The 95%% confidence interval is [%g%%, %g%%] \n', e_95L, e_95U)
% e_95L   = repmat(mean(e) - 2*std(e),[1 2]);
% e_95U   = repmat(mean(e) + 2*std(e),[1 2]);
% p1      = plot(e_95L, ax.YLim,'b--','linewidth',1);
% p2      = plot(e_95U, ax.YLim,'b--','linewidth',1);
% 
% matlab2tikz('samplingdist2.tex','showInfo', false)
% 
% title([num2str(n),'个样本时relative bias的分布']);