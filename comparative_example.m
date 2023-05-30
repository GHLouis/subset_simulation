clear; close all; clc;
figpath = 'figures\';
probname = 'comparative example';
% xlabelname = '$Y$ - Moment(kNm)';
xlabelname = '$Y$ - System response';
%% # problem satete
%{
$Y=\max \{X_1,X_2,...,X_{3000}\}$, where $X_i$s follow normal dist with mean
of 40 and std of 5.
we want to evaluate the failure probability of $Y$.
%}

% define parameters of the random variables
d      = 1;          % number of dimensions
% pi_pdf = makedist('gev','k',-0.1,'mu',3500,'sigma',265);
% pi_pdf = makedist('gev','k',0.1,'mu',3500,'sigma',265);
% pi_pdf = makedist('gev','k',0,'mu',3500,'sigma',265);
% pi_pdf = makedist('norm','mu',0,'sigma',1);
% pi_pdf = makedist('exp','mu',10);

% limit state function
beta = pi_pdf(1).icdf(exp(log(1-1/75)/(250)));
pf_ex = 1 - exp(log(1-1/75)/(250));

g = @(u) -u2x(pi_pdf,u) + beta;

% define parameters of Subset simulation
N  = 500;   % Number of samples for each level
p0 = 0.1;   % Probability of each subset, chosen adaptively
Np = 20;    % No. of repetiotion
% start subset simulation
analysis_ss
%% post process
% % opc = 'a'; % CDF of LSF
% % opc = 'b'; % CDF of response
% opc = 'c'; % CCDF of response
% % opc = 'none'; % none
% post_process_ss
post_process_ss_cv
% ylim([60 80])
%% # Start simulation
x = random(pi_pdf,[Na, Np]);
% fprintf('# *** Multiple runs ***\n')
q_hat = zeros(Np,1);
beta_hat = zeros(Np,1);
for i = 1:Np
   pd_fit        = fitdist(x(:,i), pi_pdf.DistributionName);
%    pd_fit        = ERADist('GEV','DATA',x(:,i));
   q_hat(i)      = icdf(pd_fit, 1-pf_ex);
   beta_hat(i) 	 = norminv(cdf(pd_fit,beta));
end
%%
% 特征值对比
error_ci_plot([b_SuS q_hat],beta,{'SubSim',pi_pdf.DistributionName})
xlabel('Methods','Interpreter','Latex');
ylabel('Characteristic Value, $u$','Interpreter','Latex');
c_name = 'characteristic value versus ';
% ylim([3.5 4.3])
saveas(gcf,[figpath,p_name,c_name,pi_pdf.DistributionName,'.png'])
savefig(gcf,[figpath, p_name, c_name, pi_pdf.DistributionName],'compact')
matlab2tikz([figpath,p_name,c_name,pi_pdf.DistributionName,'.tex'],'showInfo',...
    false,'checkForUpdates',false);
% 可靠指标对比
error_ci_plot([beta_SuS beta_hat],beta_ex,{'SubSim',pi_pdf.DistributionName})
xlabel('Methods','Interpreter','Latex');
ylabel('$\Phi^{-1}(P_f)$','Interpreter','Latex');
c_name = 'reliability index versus ';
% ylim([3.5 4.3])
saveas(gcf,[figpath,p_name,c_name,pi_pdf.DistributionName,'.png'])
savefig(gcf,[figpath, p_name, c_name, pi_pdf.DistributionName],'compact')
matlab2tikz([figpath,p_name,c_name,pi_pdf.DistributionName,'.tex'],'showInfo',...
    false,'checkForUpdates',false);
%% Nested function
% It is suggested to define the probability ouside the subset simulation function
% "SuS"
function x = u2x(distr,u)
x = distr(1).icdf(normcdf(u));
end