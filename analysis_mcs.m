M0      = zeros(N_MCS,1);
dim     = size(pi_pdf,1);
U       = zeros(N_MCS,dim);
X       = cell(N_MCS,1);
tic
for i = 1:N_MCS
    U(i,:)  = randn(1,dim);
    M0(i,:) = g(U(i,:));
end
toc
Pf_MCS = sum(M0<0)/N_MCS;
delta_MCS = sqrt((1-Pf_MCS)/(Pf_MCS*N_MCS));
fprintf('# MCS result: \n');
fprintf('\n MCS Pf: %g', Pf_MCS);
fprintf('\n MCS Pf COV: %g \n', delta_MCS);

M = -M0 + beta;
dist = optimal_commom_dist(M);

p_name = [probname,'_'];
c_name = 'block max fitting';
% % exportgraphics(gcf, [figpath, p_name, c_name, '.png'])
% % exportgraphics(gcf, [figpath, p_name, c_name, '.pdf'])
saveas(gcf,[figpath,p_name,c_name,'.png'])
savefig(gcf,[figpath, p_name, c_name],'compact')

save([probname,'_samples','.mat'],'M')