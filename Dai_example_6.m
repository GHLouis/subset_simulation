clear; close all; clc;
figpath = 'figures\';
probname = 'Dai example 6';
% xlabelname = '$Y$ - Moment(kNm)';
xlabelname = '$Y$ - System response';
%% # problem satete
%{
cdf of mixture distribution
%}
p = [0.5 0.3 0.2];
k =  [-0.05;nan;nan]; % gev,norm,norm
sigma = [80;120;105];
mu = [60; 215;190];

% define parameters of the random variables
d0     = 3500;
d      = d0 * 3;
fun = @(x) p(1)*gevcdf(x,k(1),sigma(1),mu(1)) ...
            + p(2)*normcdf(x,mu(2),sigma(2))...
            + p(3)*normcdf(x,mu(3),sigma(3));
%         
%         fun = @(x) p(1)*gevcdf(x,k(1),sigma(1),mu(1)) ...
%             + p(2)*normcdf(x,sigma(2),mu(2))...
%             + p(3)*normcdf(x,sigma(3),mu(3));

y = exp(log(1-1/1000)/(250*d0));
f = @(x) abs(fun(x) - y);
% limit state function
beta = 1e10;
beta0 = fminsearch(f,mean(mu));
pf_ex = 1 - exp(log(1-1/1000)/(250));
thre_opc    = 'b'; % threshold for evaluating reliability index, default is 'b'

pi_pdf = cell(d,1);
pi_pdf{1,1} = makedist('Multinomial','Probabilities',p);
pi_pdf{2,1} = k;
pi_pdf{3,1} = sigma;
pi_pdf{4,1} = mu;
g = @(u) -max(u2x(pi_pdf,u)) + beta;

% define parameters of Subset simulation
N  = 500; % Number of samples for each level
p0 = 0.1;    % Probability of each subset, chosen adaptively
Np = 20;      % No. of repetiotion
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
% M0 = -M + beta;
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
[n, d] = size(u);
u = reshape(u, [d/3 3]);
x = zeros(n,1);
k = distr{1,1}.icdf(normcdf(u(:,1)));
% counts = histcounts(k);
k = k==1;
x(k) = gevinv(normcdf(u(k,2)),distr{2,1}(1),distr{3,1}(1),distr{4,1}(1));
k = k==2;
x(k) = norminv(normcdf(u(k,2)),distr{4,1}(2),distr{3,1}(2));
k = k==3;
x(k) = norminv(normcdf(u(k,2)),distr{4,1}(3),distr{3,1}(3));
end