function [Pf_SuS,delta_SuS,b,Pf,b_line,Pf_line,samplesU,samplesX] = SuS(N, p0, g_fun, distr)
%% Subset Simulation function (standard Gaussian space)
%{
---------------------------------------------------------------------------
Created by:
Felipe Uribe
Engineering Risk Analysis Group
Technische Universitat Munchen
www.bgu.tum.de/era
Contact: Antonios Kamariotis (antonis.kamariotis@tum.de)
---------------------------------------------------------------------------
Current version 2020-10
* Adaptation to new ERANataf class
---------------------------------------------------------------------------
Comments:
*Express the failure probability as a product of larger conditional failure
 probabilities by introducing intermediate failure events.
*Use MCMC based on the modified Metropolis-Hastings algorithm for the
 estimation of conditional probabilities.
*p0, the prob of each subset, is chosen 'adaptively' to be in [0.1,0.3]
---------------------------------------------------------------------------
Input:
* N     : Number of samples per level
* p0    : Conditional probability of each subset
* g_fun : limit state function
* distr : Nataf distribution object or
          marginal distribution object of the input variables
* alg   : Sampling algorithm
          - 'acs' : Adaptive Conditional Sampling
          - 'mma' : Component-wise Metropolis Algorithm
---------------------------------------------------------------------------
Output:
* Pf_SuS    : failure probability estimate of subset simulation
* delta_SuS : coefficient of variation estimate of subset simulation
* b         : intermediate failure levels
* Pf        : intermediate failure probabilities
* b_line    : limit state function values
* Pf_line   : failure probabilities corresponding to b_line
* samplesU  : samples in the Gaussian standard space for each level
* samplesX  : samples in the physical/original space for each level
---------------------------------------------------------------------------
Based on:
1."Estimation of small failure probabilities in high dimentions by SubSim"
   Siu-Kui Au & James L. Beck.
   Probabilistic Engineering Mechanics 16 (2001) 263-277.
2. "MCMC algorithms for subset simulation"
   Papaioannou et al.
   Probabilistic Engineering Mechanics 41 (2015) 83-103.
---------------------------------------------------------------------------
%}
if (N*p0 ~= fix(N*p0)) || (1/p0 ~= fix(1/p0))
    error('N*p0 and 1/p0 must be positive integers. Adjust N and p0 accordingly');
end

%% transform to the standard Gaussian space
% if any(strcmp('Marginals',fieldnames(distr))) == 1   % use Nataf transform (dependence)
%     n   = length(distr.Marginals);    % number of random variables (dimension)
%     u2x = @(u) distr.U2X(u);          % from u to x
%     
%     % else
%     %     n   = length(distr);                    % number of random variables (dimension)
%     %     u2x = @(u) u;   % from u to x
%     
% else   % use distribution information for the transformation (independence)
%     % Here we are assuming that all the parameters have the same distribution !!!
%     % Adjust accordingly otherwise
%     n   = length(distr);                    % number of random variables (dimension)
%     %     u2x = @(u) distr(1).icdf(normcdf(u));   % from u to x
%     u2x = @(u) u;   % from u to x
% end
n   = length(distr);                    % number of random variables (dimension)
u2x = @(u) u;   % from u to x
%% LSF in standard space
G_LSF = @(u) g_fun(u2x(u));

%% Initialization of variables and storage
j      = 1;                % initial conditional level
Nc     = N*p0;             % number of markov chains
Ns     = 1/p0;             % number of samples simulated from each Markov chain
lambda = 0.6;              % recommended initial value for lambda
max_it = 6;                % maximum number of iterations
%
geval = zeros(1,N);        % space for the LSF evaluations
gsort = zeros(max_it,N);   % space for the sorted LSF evaluations
delta = zeros(max_it,1);   % space for the coefficient of variation
nF    = zeros(max_it,1);   % space for the number of failure point per level
prob  = zeros(max_it,1);   % space for the failure probability at each level
b     = zeros(max_it,1);   % space for the intermediate leveles

%% SuS procedure
% initial MCS stage
u_j = randn(n,N)';     % initial samples in the standard space
% fprintf('Evaluating LSF function:\t');
for i = 1:N
    geval(i) = G_LSF(u_j(i,:));
%     
%     Y = -geval + 20450;
end
% fprintf('OK! \n');

% SuS stage
while true
    % sort values in ascending order
    [gsort(j,:), idx] = sort(geval);
    
    % order the samples according to idx
    u_j_sort          = u_j(idx,:);   % order the samples
    samplesU.total{j} = u_j_sort;     % store the ordered samples
    
    % intermediate level
    b(j) = prctile(geval, p0*100);
    
    % number of failure points in the next level
    nF(j) = sum(geval <= max(b(j), 0));
    
    % assign conditional probability to the level
    if b(j) <= 0
        b(j)    = 0;
        prob(j) = nF(j)/N;
    else
        prob(j) = p0;
    end
%     fprintf('\n-Threshold intermediate level %g = %g \n', j-1, b(j));
    
    % compute coefficient of variation
    if j == 1
        delta(j) = sqrt(((1-p0)/(N*p0)));   % cov for p(1): MCS (Ref. 2 Eq. 8)
    else
        I_Fj     = reshape(geval <= b(j),Ns,Nc);         % indicator function for the failure samples
        p_j      = (1/N)*sum(I_Fj(:));                   % ~=p0, sample conditional probability
        gamma    = corr_factor(I_Fj,p_j,Ns,Nc);          % corr factor (Ref. 2 Eq. 10)
        delta(j) = sqrt( ((1-p_j)/(N*p_j))*(1+gamma) );  % coeff of variation(Ref. 2 Eq. 9)
    end
    
    % select seeds
    samplesU.seeds{j} = u_j_sort(1:nF(j),:);   % store ordered level seeds
    
    % randomize the ordering of the samples (to avoid bias)
    idx_rnd   = randperm(nF(j));
    rnd_seeds = samplesU.seeds{j}(idx_rnd,:);   % non-ordered seeds
    
    % MCMC sampling
    [u_j, geval, lambda, sigma, acc] = aCS(N, lambda, b(j), rnd_seeds, G_LSF);
%     fprintf('\t*aCS lambda = %g \t*aCS sigma = %g \t *aCS accrate = %g\n', lambda, sigma(1), acc);
    
    % next level
    j = j+1;
    
    if b(j-1) <= 0 || j-1 == max_it
        break;
    end
end
m = j-1;
samplesU.total{j} = u_j;   % store final failure samples (non-ordered)

% delete unnecesary data
if m < max_it
    gsort(m+1:end,:) = [];
    prob(m+1:end)    = [];
    b(m+1:end)       = [];
    delta(m+1:end)   = [];
end

%% probability of failure
% failure probability estimate
Pf_SuS = prod(prob);   % or p0^(m-1)*(Nf(m)/N);

% coeficient of variation estimate
delta_SuS = sqrt(sum(delta.^2));   % (Ref. 2 Eq. 12)

%% Pf evolution
Pf           = zeros(m,1);
Pf(1)        = p0;
Pf_line(1,:) = linspace(p0,1,Nc);
b_line(1,:)  = prctile(gsort(1,:),Pf_line(1,:)*100);
for i = 2:m
    Pf(i)        = Pf(i-1)*p0;
    Pf_line(i,:) = Pf_line(i-1,:)*p0;
    b_line(i,:)  = prctile(gsort(i,:),Pf_line(1,:)*100);
end
Pf_line = sort(Pf_line(:));
b_line  = sort(b_line(:));

%% transform the samples to the physical/original space
samplesX = cell(j,1);
for i = 1:j
    samplesX{i} = u2x(samplesU.total{i});
end

return;
%%END