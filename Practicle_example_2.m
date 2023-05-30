clear; close all; clc;
figpath = 'figures\';
probname = 'Practicle example 2';
% xlabelname = '$Y$ - Moment(kNm)';
xlabelname = '$Y$ - System response';
load('vehicle-load-model_15m_shearforce.mat')
%% # problem satete
%{
vehicle load effect problem, discret vehicle inflow
%}

% define parameters of the random variables
d      = 360;                   % number of dimensions, not ture here, the ture is 360;
pi_pdf = [
    repmat(makedist('Binomial',"N",1,"p",0.6), d, 1);...            % infolw, unit 0 or 1
    repmat(makedist('uniform',"lower",1,"upper",2), d, 1);...     % speed, unit m/s
    repmat(makedist('uniform',"lower",10,"upper",30), d, 1)];       % weight, unit kN

% limit state function
beta    = 3.6e4; % target threshold
g       = @(u) -fcn(u2x(pi_pdf,u),info) + beta; % limit state function
pf_ex   = 1 - exp(log(1-1/1000)/(250*24));
% pf_ex = 1 - exp(log(1-1/1000)/(365*24));
thre_opc    = 'a'; % threshold for evaluating reliability index, default is 'b'

% define parameters of Subset simulation
N  = 500;   % Number of samples for each level
p0 = 0.1;   % Probability of each subset, chosen adaptively
Np = 20;    % No. of repetiotion
% start subset simulation
analysis_ss
%% # Start simulation
% % MCS and autovariance
% N_MCS   = 1e2;
% M0      = zeros(N_MCS,3600);
% M       = zeros(N_MCS,1);
% M1      = zeros(N_MCS,1);
% dim     = size(pi_pdf,1);
% U       = zeros(N_MCS,dim);
% X       = cell(N_MCS,1);
% tic
% for i = 1:N_MCS
%     U(i,:)  = randn(1,dim);
%     X{i,1}  = u2x(pi_pdf,U(i,:));
%     % M0(i,:) = LE_autocorr(X{i,1},info);
% %     M1(i) = maxLE(X{i,1},info);
%     M(i) = maxLE_new(X{i,1},info);
% end
% toc
% % dist = optimal_commom_dist(M);
% % p_name = [probname,'_'];
% % c_name = 'block max fitting';
% % % % exportgraphics(gcf, [figpath, p_name, c_name, '.png'])
% % % % exportgraphics(gcf, [figpath, p_name, c_name, '.pdf'])
% % saveas(gcf,[figpath,p_name,c_name,'.png'])
% % savefig(gcf,[figpath, p_name, c_name],'compact')
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

% 子函数1：根据影响线计算荷载效应
function x = u2x(distr,u)
n   = length(distr)/3;
u   = reshape(u,n,3);

x = [distr(1).icdf(normcdf(u(:,1))) ...
    distr(n+1).icdf(normcdf(u(:,2))) ...
    distr(2*n+1).icdf(normcdf(u(:,3)))];

% x = [distr(1).random(n,1) ...
%     distr(n+1).icdf(normcdf(u(:,2))) ...
%     distr(2*n+1).icdf(normcdf(u(:,3)))];
end

function M = fcn(x,info)
M = maxLE_new(x,info);
end
% 子函数2_1：根据影响线计算荷载效应
% 按总重加载
function M = maxLE(vehicle,info)
% INPUT====================================================================
% vehicle:  [infolw, speed, weight], n by 3 matrix
% info:     vehicel load model, strcture body
% OUTPUT===================================================================
% M:        target load effect at each time step, 10n by 1 vector
% =========================================================================

n           =   size(vehicle,1);    % num of vehicles
t_sp        =   3600;               % sampling period
t_ar        =   0:10:(t_sp-10);     % vehicle arrival time
t_ar        =   t_ar(:);
v           =   vehicle(:,2);       % vehicle speed
P           =   vehicle(:,3);       % vehicle weight
dt          =   1;                  % time step for calculation
m           =   zeros(n,t_sp);      % load effect of vehicle i, at time j
x_vehicle   =   zeros(n,t_sp);      % position of vehicle i, at time j

% obtain load effect vehicle by vehicle
for i = 1 : n
    if vehicle(i,1)>0
        for j = (floor(t_ar(i))+dt):t_sp
            % update each potential car's location
            x_vehicle(i,j) = x_vehicle(i,j) + v(i)*(j-t_ar(i));
            if and(x_vehicle(i,j) >= 0, x_vehicle(i,j) < info.L)
                % update moment contribution based on potential car's new location
                %             m(i,j) = vehicle(i,4)*info.IL(x_vehicle(i,j));
                m(i,j) = P(i)*info.Inline(floor(x_vehicle(i,j)/info.dx)+1);
            else
                % ignore this potential car for the rest of the time
                break
            end
        end
    end
end
M = max(sum(m,1));
end

% 子函数2_1：根据影响线计算荷载效应
% 按总重加载
function M = LE_autocorr(vehicle,info)
% INPUT====================================================================
% vehicle:  [infolw, speed, weight], n by 3 matrix
% info:     vehicel load model, strcture body
% OUTPUT===================================================================
% M:        target load effect at each time step, 10n by 1 vector
% =========================================================================

n           =   size(vehicle,1);    % num of vehicles
t_sp        =   3600;               % sampling period
t_ar        =   0:10:(t_sp-10);     % vehicle arrival time
t_ar        =   t_ar(:);
v           =   vehicle(:,2);       % vehicle speed
P           =   vehicle(:,3);       % vehicle weight
dt          =   1;                  % time step for calculation
m           =   zeros(n,t_sp);      % load effect of vehicle i, at time j
x_vehicle   =   zeros(n,t_sp);      % position of vehicle i, at time j

% obtain load effect vehicle by vehicle
for i = 1 : n
    if vehicle(i,1)>0
        for j = (floor(t_ar(i))+dt):t_sp
            % update each potential car's location
            x_vehicle(i,j) = x_vehicle(i,j) + v(i)*(j-t_ar(i));
            if and(x_vehicle(i,j) >= 0, x_vehicle(i,j) < info.L)
                % update moment contribution based on potential car's new location
                %             m(i,j) = vehicle(i,4)*info.IL(x_vehicle(i,j));
                m(i,j) = P(i)*info.Inline(floor(x_vehicle(i,j)/info.dx)+1);
            else
                % ignore this potential car for the rest of the time
                break
            end
        end
    end
end
M = sum(m,1);     % load effect at different step
end

function M = maxLE_new(vehicle,info)
% INPUT====================================================================
% vehicle:  [infolw, speed, weight], n by 3 matrix
% info:     vehicel load model, strcture body
% OUTPUT===================================================================
% M:        target load effect at each time step, 10n by 1 vector
% =========================================================================

n           =   size(vehicle,1);    % num of vehicles
t_sp        =   3600;               % sampling period
t_ar        =   0:10:(t_sp-10);     % vehicle arrival time
t_ar        =   t_ar(:);
v           =   vehicle(:,2);       % vehicle speed
P           =   vehicle(:,3);       % vehicle weight
dt          =   1;                  % time step for calculation

M           =   zeros(3600,1);
m           =   zeros(3600,1);      % load effect of vehicle i, at time j
x_vehicle   =   zeros(3600,1);      % position of vehicle i, at time j
        
% obtain load effect vehicle by vehicle
for i = 1 : n
    % if there is a vehicle?
    if vehicle(i,1)>0
        for j = (floor(t_ar(i))+dt):t_sp
            % update each potential car's location
            x_vehicle(j) = v(i)*(j-t_ar(i));
            if x_vehicle(j) < info.L
                % update moment contribution based on potential car's new location
                m(j) = P(i)*info.Inline(floor(x_vehicle(j)/info.dx)+1);
            else
                % ignore this potential car for the rest of the time
                break
            end
        end
        M = M + m;
        m((floor(t_ar(i))+dt):j)=0;
        x_vehicle((floor(t_ar(i))+dt):j)=0;
        
    end
end
M = max(M);
end