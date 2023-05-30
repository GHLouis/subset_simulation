clear; close all; clc;
figpath = 'figures\';
probname = 'obrien example 1';
% xlabelname = '$Y$ - Moment(kNm)';
xlabelname = '$Y$ - System response';
%% # problem satete
%{
$Y=\max \{X_1,X_2,...,X_{3000}\}$, where $X_i$s follow normal dist with mean
of 40 and std of 5.
we want to evaluate the failure probability of $Y$.
%}

% define parameters of the random variables
d      = 3000;          % number of dimensions
pi_pdf = repmat(ERADist('normal','PAR',[40,5]), d, 1);

% limit state function
beta = 1e10;
beta0 = pi_pdf(1).icdf(exp(log(1-1/1000)/(250*d)));
pf_ex = 1 - exp(log(1-1/1000)/(250));
thre_opc    = 'b'; % threshold for evaluating reliability index, default is 'b'
% beta = pi_pdf(1).icdf(exp(log(1-0.05)/(365*d*100)));
% pf_ex = 1 - exp(log(1-0.05)/(365*100));
g = @(u) -max(u2x(pi_pdf,u)) + beta;

% define parameters of Subset simulation
N  = 500;   % Number of samples for each level
p0 = 0.1;   % Probability of each subset, chosen adaptively
Np = 20;    % No. of repetiotion
% start subset simulation
analysis_ss
%% # Start simulation
% % MCS
% N_MCS   = 1e3;
% M       = zeros(N_MCS,1);
% tic
% for i = 1:N_MCS
%     u       = randn(1,d);
%     M(i)    = g(u);
% end
% toc
% Pf_MCS = sum(M<0)/N_MCS;
% delta_MCS = sqrt((1-Pf_MCS)/(Pf_MCS*N_MCS));
%% post process
% % opc = 'a'; % CDF of LSF
% % opc = 'b'; % CDF of response
% opc = 'c'; % CCDF of response
% % opc = 'none'; % none
% post_process_ss
post_process_ss_cv
% ylim([60 80])
%% Nested function
% It is suggested to define the probability ouside the subset simulation function
% "SuS"
function x = u2x(distr,u)
x = distr(1).icdf(normcdf(u));
end