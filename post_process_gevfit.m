% #2: GEV fitting...
%% modified on 20230510
% average for multiple runs
lv      = ceil(-log10(pf_ex));
pf_gev  = power(10,-(1:lv));
pf_gev  = pf_gev(:);

x_gev0 = zeros(lv,Np);
for j = 1:Np
    x = pf_sus{j};
    y = beta - b_sus{j};
    [~,ia,~] = unique(x);
    x = x(ia);
    y = y(ia);
    fun = fit(x,y,'linearinterp');
    x_gev0(:,j) = fun(pf_gev);
end
x_gev = mean(x_gev0,2);


% fit tail by GEV
fprintf('## figure 3: GEV fitting for the intermediate thresholds \n');
% gev_para = [kesi,sig,mu];
[gev_para,fval] = gev_cdf_fitting(x_gev, 1-pf_gev);
r2   = fval(2);
xlabel(xlabelname,'Interpreter','Latex');
ylabel('$-\ln(-\ln(\Pr))$','Interpreter','Latex');
legend({'Daily maximum Data','GEV'},'location','northwest');
% legend off
ax = gca;
text(ax.XLim(1) + 0.2*diff(ax.XLim),ax.YLim(1) + 0.8*diff(ax.YLim),['$R^2 =',num2str(r2,4),'$'],...
    'Interpreter','Latex');
% gevinv(1-pf_ex,gev_para(1),gev_para(2),gev_para(3))

% find mean and std of lognormal dist 
R_cov = 0.1414;
fun = @(m) findlognmean(m,R_cov,gev_para);
R_mu = fminsearch(fun,beta0);
R_std = R_cov*R_mu;

gev_para_mr = gev_para;
%% one by one fit
beta_SuS = zeros(Np,1);
R = zeros(Np,1);
for j = 1:Np
    % gev_para = [kesi,sig,mu];
    [gev_para,fval] = gev_cdf_fitting(x_gev0(:,j), 1-pf_gev, 'display', false);
    r2   = fval(2);
    R(j) = r2;
    xlabel(xlabelname,'Interpreter','Latex');
    ylabel('$-\ln(-\ln(\Pr))$','Interpreter','Latex');
    beta_FORM
    beta_SuS(j) = beta_fmc;
end
%% old ver
% % one by one fit
% lv  = ceil(-log10(pf_ex));
% ypt = power(10,-(1:lv));
% ypt = ypt(:);
% beta_SuS = zeros(Np,1);
% R = zeros(Np,1);
% for j = 1:Np
%     x = beta - b_sus{j};
%     y = pf_sus{j};
%     [~,ia,~] = unique(y);
%     x = x(ia);
%     y = y(ia);
%     fun = fit(y,x,'linearinterp');
%     xpt = fun(ypt);
% 
%     % gev_para = [kesi,sig,mu];
%     [gev_para,fval] = gev_cdf_fitting(xpt, 1-ypt, 'display', false);
%     r2   = fval(2);
%     R(j) = r2;
%     xlabel(xlabelname,'Interpreter','Latex');
%     ylabel('$-\ln(-\ln(\Pr))$','Interpreter','Latex');
%     %     legend({'Daily maximum Data','GEV'},'location','northwest');
%     %     ax = gca;
%     %     text(ax.XLim(1) + 0.08*diff(ax.XLim),mean(ax.YLim),['$R^2 =',num2str(r2,4),'$'],...
%     %         'Interpreter','Latex');
% 
%     % modified on 20230509
%     %     pf   = 1 - gevcdf(beta0,gev_para(1),gev_para(2),gev_para(3));
%     % pf   = 1 - gevcdf(b_SuS(j),gev_para(1),gev_para(2),gev_para(3));
%     % beta_SuS(j) = -norminv(pf);
%     beta_FORM
%     beta_SuS(j) = beta_fmc;
% end
% 
% 
% % average for multiple runs
% x   = beta - cell2mat(b_sus);
% pt  = linspace(min(x),max(x)); pt  = pt(:);
% pf_new0 = zeros(100,Np);
% for j = 1:Np
%     x = beta - b_sus{j};
%     y = pf_sus{j};
%     [~,ia,~] = unique(x);
%     x = x(ia);
%     y = y(ia);
%     fun = fit(x,y,'linearinterp');
%     pf_new0(:,j) = fun(pt);
% end
% pf_new = mean(pf_new0,2);
% 
% Pf_SuS_mean = mean(Pf_SuS);
% lv      = floor(-log10(Pf_SuS_mean));
% pf_gev  = [power(10,-(1:lv)) Pf_SuS_mean];
% pf_gev  = pf_gev(:);
% fun     = fit(pf_new,pt,'linearinterp');
% pt_gev  = fun(pf_gev);
% 
% % fit tail by GEV
% fprintf('## figure 3: GEV fitting for the intermediate thresholds \n');
% [gev_para,fval] = gev_cdf_fitting(pt_gev, 1-pf_gev);
% r2   = fval(2);
% xlabel(xlabelname,'Interpreter','Latex');
% ylabel('$-\ln(-\ln(\Pr))$','Interpreter','Latex');
% legend({'Daily maximum Data','GEV'},'location','northwest');
% % legend off
% ax = gca;
% text(ax.XLim(1) + 0.08*diff(ax.XLim),ax.YLim(1) + 0.8*diff(ax.YLim),['$R^2 =',num2str(r2,4),'$'],...
%     'Interpreter','Latex');