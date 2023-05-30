%% # Plot failure probability: SuS
fprintf('# Simulation result: figure plot \n');
fprintf('## figure 1: failure probability plot of %g runs \n',Np);
figure('Name','Failure probability','NumberTitle','off');
hold on
for i=1:Np
    semilogy(beta - b_sus{i}, pf_sus{i})
end

figure_size = [6,3,3,2.25]; adjfig

set(gca,'yscale','log')
set(gca,'box','on')
ax = gca; % current axes
semilogy(ax.XLim,[pf_ex pf_ex],'k--')
xlabel(xlabelname,'Interpreter','Latex');
ylabel('Failure probability, $P_f$','Interpreter','Latex');


p_name = [probname,'_'];
c_name='failure prob plot';

saveas(gcf,[figpath,p_name,c_name,'.png'])
savefig(gcf,[figpath, p_name, c_name],'compact')
matlab2tikz([figpath,p_name,c_name,'.tex'],'showInfo', false,'checkForUpdates',false,'standalone',true)        
%% estimation characteristic value and reliability
b_SuS = zeros(Np,1);
for j = 1:Np
    x = beta - b_sus{j};
    y = pf_sus{j};
    [~,ia,~] = unique(y);
    x = x(ia);
    y = y(ia);
    fun = fit(y,x,'linearinterp');
    b_SuS(j) = fun(pf_ex);
end

switch thre_opc
    case 'a'
        beta0 = mean(b_SuS);
        b_ex = nan;
    case 'b'
        % beta0 = beta;
        b_ex = beta0;
end

% #1: interpolation
% beta_SuS = zeros(Np,1);
% for j = 1:Np
%     x = beta - b_sus{j};
%     y = pf_sus{j};
%     [~,ia,~] = unique(x);
%     x = x(ia);
%     y = y(ia);
%     fun = fit(x,y,'linearinterp');
%     beta_SuS(j) = -norminv(fun(beta0));
% end

% #2: GEV fitting...
post_process_gevfit
figure_size = [6,3,3,2.25]; adjfig
p_name = [probname,'_'];
c_name='gevfit';
matlab2tikz([figpath,p_name,c_name,'.tex'],'showInfo', false,'checkForUpdates',false,'standalone',true)

% exportgraphics(gcf,[figpath,p_name,c_name,'.png'])
% exportgraphics(gcf,[figpath,p_name,c_name,'.pdf'])
saveas(gcf,[figpath,p_name,c_name,'.png'])
savefig(gcf,[figpath, p_name, c_name],'compact')

% #3: Gumbel fitting...
% post_process_gevfit
%% Show p_f results and write into texfile
b_SuS_mean = mean(b_SuS);
b_SuS_cov = std(b_SuS)/mean(b_SuS);
b_er_100= (b_SuS_mean - b_ex)/b_ex*100;

fprintf('# Simulation result: numerical value \n');
fprintf(' Chracteristic value \n');
disp(table(b_ex,b_SuS_mean,b_SuS_cov,b_er_100));

beta_SuS_mean     = mean(beta_SuS);
beta_SuS_cov      = std(beta_SuS)/mean(beta_SuS);
beta_ex = -norminv(1e-6);
beta_er_100  = (beta_SuS_mean - beta_ex)/beta_ex*100;

fprintf(' Reilability index \n');
disp(table(beta_ex,beta_SuS_mean,beta_SuS_cov,beta_er_100));

N1 = length(b);
N2 = length(cell2mat(b));
Na = round((N1*N + (N2-N1)*(1-p0)*N)/Np);
fprintf('***Average NO. samples: %d ***\n\n', Na);

fid = fopen(['../report\',probname,'.tex'], 'w','n','UTF-8');
% fprintf(fid, '\\documentclass{article}\n');
% fprintf(fid, '\\begin{document}\n');
fprintf(fid, 'The mean of the characteristic value is %.2f, the relative bias is %.2f\\%%, ',...
    b_SuS_mean,b_er_100);
fprintf(fid, 'the mean of the reliability index is %.2f, the relative bias is %.2f\\%%. ',...
    beta_SuS_mean,beta_er_100);
fprintf(fid, 'The average number of samples over the %d runs is %d.\n\n'...
    , Np, Na);

fprintf(fid, '\\begin{table}[tbph]\n');
fprintf(fid, ['\\caption{Simulation result of ',probname,' (%d runs)}\n'],Np);
fprintf(fid, '\\begin{centering}\n');
fprintf(fid, '\\begin{tabular}{ccccc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Item & Exact value & Mean & COV & Relative bias (\\%%)\\tabularnewline\n');
fprintf(fid, '\\midrule\n');
fprintf(fid, 'Characteristic value & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',b_ex, b_SuS_mean,b_SuS_cov,b_er_100);
fprintf(fid, 'Reliability index & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',beta_ex, beta_SuS_mean,beta_SuS_cov,beta_er_100);
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\par\\end{centering}\n');
fprintf(fid, '\\end{table}\n\n');
%% plot figures
% figure
% hold on
% boxplot(b_SuS,'Whisker',10)
error_ci_plot(b_SuS,b_ex)
figure_size = [6,3,4,3]; adjfig
xlabel('SubSim','Interpreter','Latex');
ylabel('Characteristic Value, $u$','Interpreter','Latex');


% ylim('auto')
p_name = [probname,'_'];
c_name='characteristic value boxplot';
matlab2tikz([figpath,p_name,c_name,'.tex'],'showInfo', false,'checkForUpdates',false,'standalone',true)

fprintf(fid, '\\begin{figure}[tbph]\n');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, ['\\input{',['../code/figures/',p_name,c_name],'}\n']);
fprintf(fid, '\\end{center}\n');
fprintf(fid, ['\\caption{','Characteristic value estimation of ',probname,'}\n']);
fprintf(fid, '\\end{figure}\n\n');


% exportgraphics(gcf,[figpath,p_name,c_name,'.png'])
% exportgraphics(gcf,[figpath,p_name,c_name,'.pdf'])
saveas(gcf,[figpath,p_name,c_name,'.png'])
savefig(gcf,[figpath, p_name, c_name],'compact')

% figure
% hold on
% boxplot(beta_SuS,'Whisker',10)
error_ci_plot(beta_SuS,beta_ex)
figure_size = [6,3,4,3]; adjfig
xlabel('SubSim','Interpreter','Latex');
ylabel('$\Phi^{-1}(P_f)$','Interpreter','Latex');
% ylim('auto')
p_name = [probname,'_'];
c_name ='reliability index boxplot';
matlab2tikz([figpath,p_name,c_name,'.tex'],'showInfo', false,'checkForUpdates',false,'standalone',true)

fprintf(fid, '\\begin{figure}[tbph]\n');
fprintf(fid, '\\begin{center}\n');
fprintf(fid, ['\\input{',['../code/figures/',p_name,c_name],'}\n']);
fprintf(fid, '\\end{center}\n');
fprintf(fid, ['\\caption{','Reliability index estimation of ',probname,'}\n']);
fprintf(fid, '\\end{figure}\n\n');

% exportgraphics(gcf,[figpath,p_name,c_name,'.png'])
% exportgraphics(gcf,[figpath,p_name,c_name,'.pdf'])
saveas(gcf,[figpath,p_name,c_name,'.png'])
savefig(gcf,[figpath, p_name, c_name],'compact')

% % normal fit check
% figure
% h = histfit(b_SuS,floor(sqrt(Np)),'normal');
% h(1).FaceColor = 'none';
% xlabel('Estimation')
% ylabel('PDF')
% 
% p_name = [probname,'_'];
% c_name = 'normal-fitting-check';
% % exportgraphics(gcf, [figpath, p_name, c_name, '.png'])
% % exportgraphics(gcf, [figpath, p_name, c_name, '.pdf'])
% saveas(gcf,[figpath,p_name,c_name,'.png'])
% savefig(gcf,[figpath, p_name, c_name],'compact')
% 
% dist = optimal_commom_dist(b_SuS);

% fprintf(fid, '\\end{document}\n');
fclose(fid);