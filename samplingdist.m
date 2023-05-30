% 直接抽样法
clear;clc;close all
pd = makedist('normal','mu',0,'sigma',1);
q = 0.95;
% confidence interval
alpha = 1 - 0.95;
k = abs(norminv(alpha/2));
% 误差上限
EL = 0.01;
fprintf('具有%g%%保证率下相对误差不超过%g%%的分位值估计 \n', 100*(1-alpha), 100*EL)
%% 理论值
Q_mean = icdf(pd,q);
N = (k/EL)^2 * (sqrt(q*(1-q))/(pdf(pd, Q_mean) * Q_mean))^2;
N = ceil(N);
Q_std = sqrt(q*(1-q)/N)/normpdf(norminv(q));
Q_cov = Q_std/Q_mean;
fprintf('# 理论计算 \n')
fprintf('分位值为 %g \n', Q_mean)
% fprintf('理论上变异系数为 %g \n', Q_cov)
fprintf('需要的总样本数为 %g \n', N)
%% 基于Bootstrap方法估计分位值的标准差
% Bootstrap方法估计分位值
n       = 1000;
Q_std   = sqrt(q*(1-q)/n)/normpdf(norminv(q));
Q_cov   = Q_std/Q_mean;

x = random(pd,[n, 1]);
B = 1000; % bootstrap样本数量
q_hat = zeros(B, 1); % 存储所有分位数的向量
for i = 1:B
  boot_sample   = datasample(x, n, 'Replace', true);
  q_hat(i)      = quantile(boot_sample, q);
end
% 估计值的变异系数
c = std(q_hat)/mean(q_hat);
% 根据均值抽样分布，确定需要重复的次数
N = ceil((k/EL*c)^2);
fprintf('# Bootstrap \n')
fprintf('给定%g个样本 \n', n)
fprintf('理论上变异系数为%g \n',Q_cov)
fprintf('Bootstrap方法估计的变异系数为 %g \n', c)
fprintf('因此需要重复%g次 \n', N)
fprintf('需要的总样本数为 %g \n', n*N)
%% result
x       = random(pd,[n, N]);
q_hat   = quantile(x(:,1), q);
fprintf('# Single run \n')
Q_es_mean   = mean(q_hat);
Q_es_std    = std(q_hat);
Q_es_cov    = std(q_hat)/mean(q_hat);
er_100      = (Q_es_mean - Q_mean)/Q_mean*100;
fprintf('Estimated quantile is %g\n', Q_es_mean)
fprintf('Relative bias %g%% \n', er_100)

fprintf('# *** Multiple runs ***\n')
q_hat   = quantile(x, q);
e       = (q_hat-Q_mean)/Q_mean*100;
dist    = histfit_pdf(e(:),'norm');
xlabel('Relative bias (%)')
ylabel('PDF')
e_mean  = mean(e);
e_std   = std(e);
e_cov   = e_mean/e_std;
% er_100  = (Q_es_mean - Q_mean)/Q_mean*100;
fprintf('Mean of quantile is %g \n', mean(q_hat))
fprintf('Mean of relative bias is %g%% \n', e_mean)
fprintf('COV of relative bias is %g \n', e_cov)

ax      = gca;
e_95L   = mean(e) - 2*std(e);
e_95U   = mean(e) + 2*std(e);
% fprintf('The 95%% confidence interval is [%g%%, %g%%] \n', e_95L, e_95U)
e_95L   = repmat(mean(e) - 2*std(e),[1 2]);
e_95U   = repmat(mean(e) + 2*std(e),[1 2]);
p1      = plot(e_95L, ax.YLim,'b--','linewidth',1);
p2      = plot(e_95U, ax.YLim,'b--','linewidth',1);