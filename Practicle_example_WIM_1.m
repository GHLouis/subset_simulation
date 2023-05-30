clear; close all; clc;
figpath = 'figures\';
probname = 'Practicle example WIM 1';
% xlabelname = '$Y$ - Moment(kNm)';
xlabelname = '$Y$ - System response';
load('G107_trafficflow_Dr_1_Lane_1_WIM-vehicle-load-model_35m_moment.mat')
%% # problem satete
%{
vehicle load effect problem, discret vehicle inflow
%}
d           = sum(info.q);
pi_pdf      = cell(3*d,1);
pi_pdf{1}   = info.pd_cls;              % vehicel type
pi_pdf{2}   = info.pd_d;                % inter-vehicle distance
pi_pdf{3}   = info.pdw_new;             % vehicle weight

% limit state function
beta    = 3.6e10;                                           % target threshold
g       = @(x) -obj_fcn(u2x(pi_pdf,x,info),info) + beta;    % limit state function
pf_ex   = 1 - exp(log(1-1/1000)/(250));
% pf_ex = 1 - exp(log(1-1/100)/(365));
thre_opc    = 'a'; % threshold for evaluating reliability index, default is 'b'

% define parameters of Subset simulation
N  = 500;   % Number of samples for each level
p0 = 0.1;   % Probability of each subset, chosen adaptively
Np = 20;    % No. of repetiotion
%% start subset simulation
analysis_ss
%% vehicle load visiualization
U = samplesU.total{1, 6};
X = u2x(pi_pdf,U(round(N*0.4),:),info);
[M_X, idM] = maxLE_axle_gd(X,info);
% p_name = [probname,'_'];
% c_name='most unfavorable loading';
% 
% saveas(gcf,[figpath,p_name,c_name,'.png'])
% savefig(gcf,[figpath, p_name, c_name],'compact')
% matlab2tikz([figpath,p_name,c_name,'.tex'],'showInfo', false,'checkForUpdates',false,'standalone',true)
% 
% vehicle = X;
% x0 = [0;-cumsum(vehicle(1:end-1,4))];
% x0(2:end) = x0(2:end) - cumsum(vehicle(1:end-1,2));
% y = X(:,3);
% y0 = zeros(size(y));
%% # Start simulation
% MCS and autovariance
N_MCS   = 1e2;
M0      = zeros(N_MCS,3600);
M       = zeros(N_MCS,1);
M1      = zeros(N_MCS,1);
dim     = size(pi_pdf,1);
U       = zeros(N_MCS,dim);
X       = cell(N_MCS,1);
tic
for i = 1:N_MCS
    U(i,:)  = randn(1,dim);
    X{i,1}  = u2x(pi_pdf,U(i,:),info);
    M(i) = maxLE_axle(X{i,1},info);
end
toc
fprintf('mean of LE is %g\n',mean(M))
fprintf('std of LE is %g\n',std(M))

% dist = optimal_commom_dist(M);
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

% 子函数1：变量转化
function x = u2x(distr,u,info)
% OUTPUT===================================================================
% x: vehicle load. matrix
% = [ vehicle type(1), inter-vehicle distance(2), total weight(3), total length(4), 
% axle weight(5-10), axle spacing(11-16)];
% =========================================================================
n  = length(distr)/3;
u  = reshape(u,n,3);

% 车辆类型和车间距
cls = zeros(n,1);
dis = zeros(n,1);

cls(1:info.qcs(1),:)  = distr{1, 1}{1, 1}.icdf(normcdf(u(1:info.qcs(1),1)));
dis(1:info.qcs(1),:)  = distr{2, 1}{1, 1}.icdf(normcdf(u(1:info.qcs(1),2)));
for i = 2:24
    cls(info.qcs(i-1)+1:info.qcs(i),:)  = distr{1, 1}{i, 1}.icdf(normcdf(u(info.qcs(i-1)+1:info.qcs(i),1)));
    dis(info.qcs(i-1)+1:info.qcs(i),:)  = distr{2, 1}{i, 1}.icdf(normcdf(u(info.qcs(i-1)+1:info.qcs(i),2)));
end
dis(dis>=info.L) = info.L;

% ERADist('gamma','DATA',x)
% 车间距
% dis = distr{2, 1}.icdf(normcdf(u(:,2)));



% 车长车重
W   = zeros(n,1);
for j = 1:length(distr{3, 1})
    W(cls==j,1) = distr{3, 1}(j, 1).icdf(normcdf(u(cls==j,3)));
end

x = [info.Numaxle(cls) dis W];

% axle weight and wheelbase
AS = zeros(n,6); % 轴距
AW = zeros(n,6); % 轴重
% total length
TL = zeros(n,1);
% % 1 类车
% i = 1;
% AW(cls==i,1:2)  = x(cls==i,3).*info.wratio(i,1:2);
% AS(cls==i,2:2)  = repmat(info.wheelbase(1,1),[sum(cls==i) 1]);
% 
% % 2~6 类车
% for i = 2:6
%     AW(cls==i,1:i)  = x(cls==i,3).*info.wratio(i,1:i);
%     AS(cls==i,2:i)  = repmat(info.wheelbase(i,1:i-1),[sum(cls==i) 1]);
% end

for i = 1:info.Npdw
    AW(cls==i,:)  = x(cls==i,3).*info.wratio(i,:);
    AS(cls==i,2:6)= info.one(cls==i).*info.wheelbase(i,:);
    % total length
    TL(cls==i,:) = info.VL(i,:);
end

x = [x TL AW AS];
end

function M = obj_fcn(x,info)
% 计算荷载效应
% M = maxLE(x,info); % total weight
M = maxLE_axle(x,info); % axle weight
end

% 子函数2_1：根据影响线计算荷载效应
% 按总重加载
function M = maxLE(vehicle,info)
% modified on 20230328
% INPUT=================================================
% vehicle:  vehicle type, inter-vehicle distance, weight, and etc. matrix
% info:     vehicel load model. strcture body
% OUTPUT================================================
% M:        target load effect. scalar
% ======================================================
n           =   size(vehicle,1); % num of vehicles
% dx          =   info.dx;                % step for calculation
dx          =   1; % step for calculation

% initial position of front axle of each vechile
% assume vehicle length is 8m, center of weight is at the middle of
% vehicle.
x0 = [0;-cumsum(vehicle(1:end-1,4))];
x0(2:end) = x0(2:end) - cumsum(vehicle(1:end-1,2));
n_dx = ceil((-x0(end)+info.L)/dx); % No. of vehicle moving steps
M = zeros(n_dx,1);
m = zeros(n_dx,1); % load effect of vehicle i, at time j


% obtain load effect vehicle by vehicle
for i = 1 : n
    if vehicle(i,3) > 300 % 超过30吨的车，才有必要算
        j = -ceil(x0(i)) + (dx:info.L);
        x_vehicle = x0(i) + j;
        % vehicle by vehicle
        m(j) = vehicle(i,3)*info.Inline(floor(x_vehicle/info.dx)+1); 
        M(j) = M(j) + m(j);
        m(j)=0;
    end
end
M = max(M);
end

function M = maxLE_axle(vehicle,info)
% modified on 20230503
% INPUT=================================================
% vehicle:  vehicle type, inter-vehicle distance, weight, and etc. matrix
% info:     vehicel load model. strcture body
% OUTPUT================================================
% M:        target load effect. scalar
% ======================================================
n           =   size(vehicle,1); % num of vehicles
% dx          =   info.dx;                % step for calculation
dx          =   1; % step for calculation

% initial position of front axle of each vechile
% assume vehicle length is 8m, center of weight is at the middle of
% vehicle.
x0 = [0;-cumsum(vehicle(1:end-1,4))];
x0(2:end) = x0(2:end) - cumsum(vehicle(1:end-1,2));
n_dx = ceil((-x0(end)+info.L)/dx); % No. of vehicle moving steps
M = zeros(n_dx,1);
m = zeros(n_dx,1);

% obtain load effect vehicle by vehicle
for i = 1 : n
    if vehicle(i,3) > 300 % 超过30吨的车，才有必要算
        j = -ceil(x0(i)) + (dx:info.L);
        x_vehicle = x0(i) + j;
        
%         % vehicle by vehicle
%         m(j) = vehicle(i,3)*info.Inline(floor(x_vehicle/info.dx)+1);
        
        % current vehicle, axle by axle
        x_k = x_vehicle(:) - cumsum(vehicle(i,11:11+vehicle(i,1)-1));
        x_k = min(max(x_k,0),info.L);
        m(j) = info.Inline(floor(x_k/info.dx)+1) * transpose(vehicle(i,5:5+vehicle(i,1)-1));
        
        M(j) = M(j) + m(j);
        m(j)=0;
    end
end
M = max(M);
end