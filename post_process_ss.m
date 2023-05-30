%% Show p_f results
fprintf('# Simulation result: numerical value \n');
fprintf('***SINGLE RUN***');
fprintf('\n SuS Pf: %g', Pf_SuS(1));
fprintf('\n SuS Pf COV: %g \n', delta_SuS(1));
Pf_SuS_mean     = mean(Pf_SuS);
Pf_SuS_cov      = std(Pf_SuS)/mean(Pf_SuS);

% exact solution if exist
er_100  = (Pf_SuS_mean - pf_ex)/pf_ex*100;

fprintf('***MULTIPLE RUNS ***\n');
fprintf(' Failure probability \n');
disp(table(pf_ex,Pf_SuS_mean,Pf_SuS_cov,er_100));

fprintf(' Reliability index \n');
b_ex = norminv(1-pf_ex);
b_SuS = norminv(1-Pf_SuS);
b_SuS_mean = mean(b_SuS);
b_SuS_cov = std(norminv(1-Pf_SuS))/mean(norminv(1-Pf_SuS));
er_100= (b_SuS_mean - b_ex)/b_ex*100;
disp(table(b_ex,b_SuS_mean,b_SuS_cov,er_100));

N1 = length(b);
N2 = length(cell2mat(b));
fprintf('***Average NO. samples: %g ***\n\n', (N1*N + (N2-N1)*(1-p0)*N)/Np);

% normal fit check
figure
h = histfit(Pf_SuS,floor(sqrt(Np)),'normal');
h(1).FaceColor = 'none';
xlabel('Estimation')
ylabel('PDF')

p_name = [probname,'_'];
c_name = 'normal-fitting-check';
% exportgraphics(gcf, [figpath, p_name, c_name, '.png'])
% exportgraphics(gcf, [figpath, p_name, c_name, '.pdf'])
saveas(gcf,[figpath,p_name,c_name,'.png'])
savefig(gcf,[figpath, p_name, c_name],'compact')

dist = optimal_commom_dist(Pf_SuS);
%% # Plot failure probability: SuS
switch opc
    case 'a'
        % multiple runs
        fprintf('# Simulation result: figure plot \n');
        fprintf('## figure 1: failure probability plot of %g runs \n',Np);
        figure('Name','Failure probability','NumberTitle','off');
        hold on
        for i=1:Np
            semilogy(b_sus{i}, pf_sus{i})
        end
        set(gca,'yscale','log')
        set(gca,'box','on')
        ax = gca; % current axes
        semilogy(ax.XLim,[pf_ex pf_ex],'k--')
        semilogy([0 0],ax.YLim,'k--')
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
        x   = cell2mat(b_sus);
        pt  = linspace(min(x),max(x)); pt  = pt(:);
        pf_new0 = zeros(100,Np);
        for j = 1:Np
            %     x = beta-b_sus{j};
            x = b_sus{j};
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
        %
        pf_new = mean(pf_new0,2);
        
        fprintf('## figure 2: failure probability plot averaged for the %g runs \n',Np);
        figure('Name','Failure probability (averaged)','NumberTitle','off');
        hold on
        semilogy(pt, pf_new,'LineWidth',1.5)
        % plot(xi,Pf_mcs,'-.','LineWidth',1.5)
        set(gca,'yscale','log')
        set(gca,'box','on')
        xlabel(xlabelname,'Interpreter','Latex');
        ylabel('Failure probability, $P_{F}(u)$','Interpreter','Latex');
        grid on
        % legend({'SubSim','MCS'})
        legend({'SubSim'})
        xlim(ax.XLim)
        ylim(ax.YLim)
    case 'b'
        % multiple runs
        fprintf('# Simulation result: figure plot \n');
        fprintf('## figure 1: failure probability plot of %g runs \n',Np);
        figure('Name','Failure probability','NumberTitle','off');
        hold on
        for i=1:Np
            semilogy(b_sus{i} + beta, pf_sus{i})
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
        x   = cell2mat(b_sus) + beta;
        pt  = linspace(min(x),max(x)); pt  = pt(:);
        pf_new0 = zeros(100,Np);
        for j = 1:Np
            x = b_sus{j} + beta;
            y = pf_sus{j};
            [~,ia,~] = unique(x);
            x = x(ia);
            y = y(ia);
            fun = fit(x,y,'linearinterp');
            pf_new0(:,j) = fun(pt);
        end
        pf_new = mean(pf_new0,2);
        
        fprintf('## figure 2: failure probability plot averaged for the %g runs \n',Np);
        figure('Name','Failure probability (averaged)','NumberTitle','off');
        hold on
        semilogy(pt, pf_new,'LineWidth',1.5)
        % plot(xi,Pf_mcs,'-.','LineWidth',1.5)
        set(gca,'yscale','log')
        set(gca,'box','on')
        xlabel(xlabelname,'Interpreter','Latex');
        ylabel('Failure probability, $P_{F}(u)$','Interpreter','Latex');
        grid on
        % legend({'SubSim','MCS'})
        legend({'SubSim'})
        xlim(ax.XLim)
        ylim(ax.YLim)
    case 'c'
        % multiple runs
        fprintf('# Simulation result: figure plot \n');
        fprintf('## figure 1: failure probability plot of %g runs \n',Np);
        figure('Name','Failure probability','NumberTitle','off');
        hold on
        for i=1:Np
            semilogy(beta - b_sus{i}, pf_sus{i})
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
        x   = beta - cell2mat(b_sus);
        pt  = linspace(min(x),max(x)); pt  = pt(:);
        pf_new0 = zeros(100,Np);
        for j = 1:Np
            x = beta - b_sus{j};
            y = pf_sus{j};
            [~,ia,~] = unique(x);
            x = x(ia);
            y = y(ia);
            fun = fit(x,y,'linearinterp');
            pf_new0(:,j) = fun(pt);
        end
        pf_new = mean(pf_new0,2);
        
        fprintf('## figure 2: failure probability plot averaged for the %g runs \n',Np);
        figure('Name','Failure probability (averaged)','NumberTitle','off');
        hold on
        semilogy(pt, pf_new,'LineWidth',1.5)
        % plot(xi,Pf_mcs,'-.','LineWidth',1.5)
        set(gca,'yscale','log')
        set(gca,'box','on')
        xlabel(xlabelname,'Interpreter','Latex');
        ylabel('Failure probability, $P_{F}(u)$','Interpreter','Latex');
        grid on
        % legend({'SubSim','MCS'})
        legend({'SubSim'})
        xlim(ax.XLim)
        ylim(ax.YLim)
        
        %% # GEV fitting-20221213 the latest version
        lv      = floor(-log10(Pf_SuS_mean));
        pf_gev  = [power(10,-(1:lv)) Pf_SuS_mean];
        pf_gev  = pf_gev(:);
        fun     = fit(pf_new,pt,'linearinterp');
        pt_gev  = fun(pf_gev);
        
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
end
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