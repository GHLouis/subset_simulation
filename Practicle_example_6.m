%% # problem satete
%{

vehicle load effect problem, continuous vehicle inflow

%}
clear; close all; clc;
figpath = 'figures\';
probname = 'Practicle example 6';
% xlabelname = '$Y$ - Moment(kNm)';
xlabelname = '$Y$ - System response';
load('vehicle-load-model_1000.mat')
% # Initialization and parameters definition
Np          = 100;              % No. of repetiotion, 100 for default
Pf_SuS      = zeros([Np 1]);    % failure probability
delta_SuS   = zeros([Np 1]);    % cov of each run
b           = cell([Np 1]);     % intermediate failure levels
Pf          = cell([Np 1]);     % intermediate failure probabilities
b_sus       = cell([Np 1]);     % limit state function values
pf_sus      = cell([Np 1]);     % failure probabilities corresponding to b_sus

% define parameters of the random variables
d      = 360;                   % number of dimensions, not ture here, the ture is 360;
pi_pdf = [
    repmat(makedist('exp',10), d, 1);...                             % infolw, unit 0 or 1
    repmat(makedist('norm',"mu",10,"sigma",2), d, 1);...            % speed, unit m/s
    repmat(makedist('uniform',"lower",10,"upper",20), d, 1)];       % weight, unit kN

% limit state function
beta = 7.5e4; % target threshold
g    = @(u) -fcn(u2x(pi_pdf,u),info) + beta; % limit state function

% define parameters of Subset simulation
N  = 4000;    % Number of samples for each level
p0 = 0.1;     % Probability of each subset, chosen adaptively
%% # Start simulation
% MCS and autovariance
N_MCS   = 1e3;
M0      = zeros(N_MCS,3600);
M       = zeros(N_MCS,1);
M1      = zeros(N_MCS,1);
dim     = size(pi_pdf,1);
U       = zeros(N_MCS,dim);
X       = cell(N_MCS,1);
tic
for i = 1:N_MCS
    U(i,:)  = randn(1,dim);
    X{i,1}  = u2x(pi_pdf,U(i,:));
    % M0(i,:) = LE_autocorr(X{i,1},info);
%     M1(i) = maxLE(X{i,1},info);
    M(i) = maxLE_new(X{i,1},info);
end
toc
dist = optimal_commom_dist(M);
% p_name = [probname,'_'];
% c_name = 'block max fitting';
% % exportgraphics(gcf, [figpath, p_name, c_name, '.png'])
% % exportgraphics(gcf, [figpath, p_name, c_name, '.pdf'])
% saveas(gcf,[figpath,p_name,c_name,'.png'])
% savefig(gcf,[figpath, p_name, c_name],'compact')
%% 相关性分析
% rho = corr(M0);
% M_1 = mean(M0);
% M_2 = std(M0);
% M_3 = moment(M0,3)./(std(M0).^3);
% M_4 = moment(M0,4)./(std(M0).^4);
% % M_5 = moment(M0,5)./(std(M0).^5);
% 
% t = 1:size(M0,2);
% figure(1)
% hold on
% plot(t, M_1,'linewidth',1.5)
% plot(t, M_2,'--','linewidth',1.5)
% plot(t, M_3,'-.','linewidth',1.5)
% plot(t, M_4,':','linewidth',1.5)
% xlabel('Time(s)')
% ylabel('System response')
% set(gca,'yscale','log')
% set(gca,'box','on')
% grid on
% legend({'mean','std','skewness','kurtosis'})
% p_name = [probname,'_'];
% c_name = 'moments';
% % exportgraphics(gcf, [figpath, p_name, c_name, '.png'])
% % exportgraphics(gcf, [figpath, p_name, c_name, '.pdf'])
% saveas(gcf,[figpath,p_name,c_name,'.png'])
% savefig(gcf,[figpath, p_name, c_name],'compact')
% 
% figure(2);
% MM = M0(:,0.5*size(M0,2)+[0 2 5 10]);
% plotmatrix(MM)
% p_name = [probname,'_'];
% c_name = 'corr matrix';
% % exportgraphics(gcf, [figpath, p_name, c_name, '.png'])
% % exportgraphics(gcf, [figpath, p_name, c_name, '.pdf'])
% saveas(gcf,[figpath,p_name,c_name,'.png'])
% savefig(gcf,[figpath, p_name, c_name],'compact')
% 
% for i = 1:4
%     dist = optimal_commom_dist(MM(:,i));
%     close gcf
%     dist.opdist
% end
% 
% 
% % plot(t - t(5*d+1), rho(:,5*d+1),'DisplayName','rho')
% x = t - t(0.5*size(M0,2)+1); x = x(:);
% y = rho(:,0.5*size(M0,2)+1); y = y(:);
% obfun = @(a, x) exp(-(x/a).^2);
% [f,gof] = fit(x, y, obfun,'StartPoint',100);
% figure(3)
% plot(f, x, y);
% xlabel('r - t (s)')
% ylabel('Correlation coefficient')
% grid on
% text(0,0.8,['$R^2 =',num2str(gof.rsquare,4),'$'],...
%     'Interpreter','Latex');
% % xlim([-200 200])
% p_name = [probname,'_'];
% c_name = 'correlation coefficient';
% % exportgraphics(gcf, [figpath, p_name, c_name, '.png'])
% % exportgraphics(gcf, [figpath, p_name, c_name, '.pdf'])
% saveas(gcf,[figpath,p_name,c_name,'.png'])
% savefig(gcf,[figpath, p_name, c_name],'compact')


% M    = max(M0,[],2);
%%
% Subset simulation
WB = waitbar(0,'0','Name','Subset simulation running...');
tic
for i=1:Np
    [Pf_SuS(i), delta_SuS(i), b{i}, Pf{i}, b_sus{i}, pf_sus{i}, samplesU, samplesX] = SuS(N,p0,g,pi_pdf);
    waitbar(i/Np,WB,sprintf(['Finish ',num2str(round(100*i/Np)),'%%']))
end
close(WB)
time=toc;

% save the result
t = datestr(datetime('now'),'yyyymmdd');
save([probname,'_',t,'.mat'],'Pf_SuS', 'delta_SuS', 'b', 'Pf', 'b_sus',...
    'pf_sus','time','beta')
%% Show p_f results
fprintf('# Simulation result: numerical value \n');
fprintf('***SINGLE RUN***');
fprintf('\n SuS Pf: %g', Pf_SuS(1));
fprintf('\n SuS Pf COV: %g \n', delta_SuS(1));
Pf_SuS_mean     = mean(Pf_SuS);
Pf_SuS_cov      = std(Pf_SuS)/mean(Pf_SuS);

% exact solution if exist
pf_ex   = nan;
er_100  = (Pf_SuS_mean - pf_ex)/pf_ex*100;

fprintf('***MULTIPLE RUNS ***\n');
fprintf(' Failure probability \n');
disp(table(pf_ex,Pf_SuS_mean,Pf_SuS_cov,er_100));

fprintf(' Reliability index \n');
b_ex = norminv(1-pf_ex);
b_SuS_mean = mean(norminv(1-Pf_SuS));
b_SuS_cov = std(norminv(1-Pf_SuS))/mean(norminv(1-Pf_SuS));
er_100= (b_SuS_mean - b_ex)/b_ex*100;
disp(table(b_ex,b_SuS_mean,b_SuS_cov,er_100));

fprintf('***Average NO. samples: %g ***\n\n', length(cell2mat(b))/Np*N);
%% # Plot failure probability: SuS
% multiple runs
fprintf('# Simulation result: figure plot \n');
fprintf('## figure 1: failure probability plot of %g runs \n',Np);
figure('Name','Failure probability','NumberTitle','off');
hold on
for i=1:Np
    semilogy(beta-b_sus{i}, pf_sus{i})
end
set(gca,'yscale','log')
set(gca,'box','on')
ax = gca; % current axes
semilogy(ax.XLim,[pf_ex pf_ex],'k--')
semilogy([beta beta],ax.YLim,'k--')
text(1.05*ax.XLim(1),power(10,1-ceil(-log10(pf_ex))),'Exact $P_{F}(u^{\#})$',...
    'Interpreter','Latex');
text(0.95*beta,1e-1,'Exact $u^{*}$','Interpreter','Latex');
xlabel(xlabelname,'Interpreter','Latex');
ylabel('Failure probability, $P_{F}(u)$','Interpreter','Latex');

grid on
p_name = [probname,'_'];
c_name = 'failure_prob_plot1';
% exportgraphics(gcf, [figpath, p_name, c_name, '.png'])
% exportgraphics(gcf, [figpath, p_name, c_name, '.pdf'])
saveas(gcf,[figpath,p_name,c_name,'.png'])
savefig(gcf,[figpath, p_name, c_name],'compact')

% average for multiple runs
x = beta-b_sus{1};
pt = linspace(min(x),beta);
pt = pt(:);
pf_new0 = zeros(100,Np);
for j = 1:Np
    x = beta-b_sus{j};
    y = pf_sus{j};
    %     figure;plot(x,y)
    %     set(gca,'yscale','log')
    [~,ia,~] = unique(x);
    x = x(ia);
    y = y(ia);
    %     figure;plot(x,y)
    %     set(gca,'yscale','log')
    fun = fit(x,y,'linearinterp');
    pf_new0(:,j) = fun(pt);
end

pf_new = mean(pf_new0,2);

% MCS
[Pf_mcs,xi] = ecdf(M,'function',"survivor");

% % Exact solution if exist
% b_new = linspace(ax.XLim(1),beta);
% b_new = b_new(:);
% Pf_exact = zeros(1e2,1);
% for i=1:100
%     B = b_new(i);
%     fun = @(x) normcdf(-x);
%     Pf_exact(i) = fun(B);
% end
% Pf_exact_51 = 0.95*Pf_exact;
% Pf_exact_52 = 1.05*Pf_exact;

fprintf('## figure 2: failure probability plot averaged for the %g runs \n',Np);
figure('Name','Failure probability (averaged)','NumberTitle','off');
hold on
plot(pt, pf_new,'LineWidth',1.5)
plot(xi,Pf_mcs,'-.','LineWidth',1.5)
set(gca,'yscale','log')
set(gca,'box','on')
xlabel(xlabelname,'Interpreter','Latex');
ylabel('Failure probability, $P_{F}(u)$','Interpreter','Latex');
grid on
legend({'SubSim','MCS'})
xlim(ax.XLim)

try
    plot(b_new,Pf_exact,'--','LineWidth',1.5)
    plot(b_new,Pf_exact_51,'--','LineWidth',1.5)
    plot(b_new,Pf_exact_52,'--','LineWidth',1.5)
    legend({'SubSim','MCS','Exact','Exact - 5%','Exact + 5%'})
catch
    fprintf('No exact solution \n');
end


c_name='failure_prob_plot2';
% exportgraphics(gcf,[figpath,p_name,c_name,'.png'])
% exportgraphics(gcf,[figpath,p_name,c_name,'.pdf'])
saveas(gcf,[figpath,p_name,c_name,'.png'])
savefig(gcf,[figpath,p_name,c_name],'compact')
%% # GEV fitting-20221213 the latest version
lv = floor(-log10(Pf_SuS_mean));
pf_gev = [power(10,-(1:lv)) Pf_SuS_mean];
pf_gev = pf_gev(:);
fun = fit(pf_new,pt,'linearinterp');
pt_gev = fun(pf_gev);

% fit tail by GEV
fprintf('## figure 3: GEV fitting for the intermediate thresholds \n');
[gev_para,fval] = gev_cdf_fitting(pt_gev, 1-pf_gev);
R2   = fval(2);
xlabel(xlabelname,'Interpreter','Latex');
ylabel('$-\ln(-\ln(P_{F}(u)))$','Interpreter','Latex');
legend('Daily maximum Data','GEV','location','northwest');
text(pt_gev(2),-log(-log(1-pf_gev(end-1))),['$R^2 =',num2str(R2,4),'$'],...
    'Interpreter','Latex');
set(gca,'Box','on')
grid on
c_name = 'gevfiting';
saveas(gcf,[figpath,p_name,c_name,'.png'])
% exportgraphics(gcf,[figpath,p_name,c_name,'.pdf'])
savefig(gcf,[figpath,p_name,c_name],'compact')

% or fit tail by Gumbel
fprintf('## figure 4: Gumbel fitting for the intermediate thresholds \n');
[fun,gof] = gumbel_fit(pt_gev, 1-pf_gev);
R2   = gof.rsquare;
xlabel(xlabelname,'Interpreter','Latex');
ylabel('$-\ln(-\ln(P_{F}(u)))$','Interpreter','Latex');
legend('Daily maximum Data','Gumbel','location','northwest');
text(pt_gev(2),-log(-log(1-pf_gev(end-1))),['$R^2 =',num2str(R2,4),'$'],...
    'Interpreter','Latex');
set(gca,'Box','on')
grid on

c_name = 'gumbelfiting';
saveas(gcf,[figpath,p_name,c_name,'.png'])
% exportgraphics(gcf,[figpath,p_name,c_name,'.pdf'])
savefig(gcf,[figpath,p_name,c_name],'compact')

% extrapolate P_F
pf   = 1 - gevcdf(beta,gev_para(1),gev_para(2),gev_para(3));
% % extrapolate P_M
% pM   = 1 - (gevcdf(beta,gev_para(1),gev_para(2),gev_para(3)))^(365*100);
% extrapolate threshold
thre = gevinv(1 - pf_ex,gev_para(1),gev_para(2),gev_para(3));

fprintf('\n # Failure probability extrapolation result \n');
fprintf('## Estimated failure probability Pf is:  %g \n',pf);
% fprintf('## Estimated failure probability PM is:  %g \n',pM);
fprintf('## Estimated characteristic value u is:  %g \n',thre);
%% GEV fitting - the old version
% pf      = zeros(Np,1);
% pM      = zeros(Np,1);
% R2      = zeros(Np,1);
% thre    = zeros(Np,1);
% for i=1:Np
%     b{i}(end) = b_sus{i}(1);
%     [gev_para,fval] = gev_cdf_fitting(beta-b{i}, 1-Pf{i});
%     xlabel(xlabelname,'Interpreter','Latex');
%     ylabel('$-\ln(-\ln(\mathrm{P}))$','Interpreter','Latex');
%     legend('Daily maximum Data','GEV fitting','location','northwest');
%     xlim([0 12]);
%     %       [gev_para,fval]=gev_cdf_fitting(beta-b{i}(1:end-1), 1-Pf{i}(1:end-1));
%     %       [gev_para,fval]=gev_cdf_fitting(beta-b_sus{i}(1:end-1), 1-pf_sus{i}(1:end-1));
%     R2(i)   = fval(2);
%     pf(i)   = 1-gevcdf(beta,gev_para(1),gev_para(2),gev_para(3)); % extrapolate P_F
%     pM(i)   = 1-(gevcdf(beta,gev_para(1),gev_para(2),gev_para(3)))^(250*75); % extrapolate P_M
%     thre(i) = gevinv(1-1/(250*75),gev_para(1),gev_para(2),gev_para(3)); % extrapolate threshold
%     tx = beta-b{1};
%     ty = -log(-log(1-Pf{i}));
%     text(tx(1),ty(end-1),['$R^2 =',num2str(R2(i),4),'$'],'Interpreter','Latex');
%     set(gca,'Box','on')
%     set(gca,'FontName','Times New Roman')
%     grid on
% end
% c_name = 'gevfiting';
% exportgraphics(gcf,[figpath,p_name,c_name,'.png'])
% exportgraphics(gcf,[figpath,p_name,c_name,'.pdf'])
% savefig(gcf,[figpath,p_name,c_name],'compact')
%%
% figure
% hold on
% plot(0:2,zeros(1,3)+beta,'k--')
% boxplot(thre,'Whisker',10)
% ylabel('Threshold, $u$','Interpreter','Latex');
% ylim('auto')
% p_name='main_1_';
% c_name='thre_boxplot';
% exportgraphics(gcf,[figpath,p_name,c_name,'.png'])
% exportgraphics(gcf,[figpath,p_name,c_name,'.pdf'])
% savefig(gcf,[figpath,p_name,c_name],'compact')
%
% figure
% hold on
% plot(0:2,zeros(1,3)+pf_ex,'k--')
% boxplot(pf,'Whisker',10)
% ylim('auto')
% ylabel('Failure probability, $P_{F}\left(u\right)$','Interpreter','Latex');
% c_name='pf_boxplot';
% exportgraphics(gcf,[figpath,p_name,c_name,'.png'])
% exportgraphics(gcf,[figpath,p_name,c_name,'.pdf'])
% savefig(gcf,[figpath,p_name,c_name],'compact')
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
M = maxLE(x,info);
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
t_sp        =   3600;                 % sampling period   
t_ar        =   cumsum(vehicle(:,1)); % vehicle arrival time
v           =   vehicle(:,2);       % vehicle speed
P           =   vehicle(:,3);       % vehicle weight
dt          =   1;                  % time step for calculation
m           =   zeros(n,t_sp);      % load effect of vehicle i, at time j
x_vehicle   =   zeros(n,t_sp);      % position of vehicle i, at time j

% obtain load effect vehicle by vehicle
for i = 1 : n
    if t_ar(i) < t_sp
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
t_sp        =   36;                 % sampling period
t_ar        =   cumsum(vehicle(:,1)); % vehicle arrival time
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
t_ar        =   cumsum(vehicle(:,1)); % vehicle arrival time
v           =   vehicle(:,2);       % vehicle speed
P           =   vehicle(:,3);       % vehicle weight
dt          =   1;                  % time step for calculation

M           =   zeros(3600,1);
m           =   zeros(3600,1);      % load effect of vehicle i, at time j
x_vehicle   =   zeros(3600,1);      % position of vehicle i, at time j

% obtain load effect vehicle by vehicle
for i = 1 : n
    % if there is a vehicle?
    if t_ar(i)< (t_sp+1)
        for j = (floor(t_ar(i))+dt):t_sp
            % update each potential car's location
            x_vehicle(j) = v(i)*(j-t_ar(i));
            if x_vehicle(j) < info.L
                % update moment contribution based on potential car's new location
                % 这里可能非整数？？？ 需要解决
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