%%
pi_pdf = [  ERADist('lognormal','MOM',[R_mu,R_std]);
            ERADist('GEV','PAR',[gev_para(1),gev_para(2),gev_para(3)])];

% correlation matrix
Rx = eye(2);   % independent case

% object with distribution information
pi_pdf = ERANataf(pi_pdf,Rx);    % if you want to include dependence

% limit state function and its gradient in the original space
g    = @(x) x(1)-x(2);

% Solve the optimization problem of the First Order Reliability Method
% FORM using MATLAB fmincon
[~, ~, beta_fmc, ~] = FORM_fmincon(g, pi_pdf);
% beta_tivr = beta_fmc;
% show p_f results
% print results
% fprintf(' Reliability index = %g --- Failure probability = %g\n\n',beta_fmc,Pf_fmc);
% fprintf(' Design points: %g, %g\n\n',x_star_fmc(1),x_star_fmc(2));
% fprintf(' Performance function: %g\n\n',g(x_star_fmc));
% %%
% x = linspace(icdf(pi_pdf.Marginals(2, 1),1e-4) ,icdf(pi_pdf.Marginals(1, 1),1-1e-4),200);
% y1 = pdf(pi_pdf.Marginals(1, 1),x); % R
% y2 = pdf(pi_pdf.Marginals(2, 1),x); % S
% figure
% hold on
% plot(x,y1)
% plot(x,y2)
% legend({'R','S'})
% xlabel(xlabelname,'Interpreter','latex')
% ylabel('PDF')
% box("on")
%% 
% function y = findlognmean(m,c,u)
% s = c*m;
% pd = ERADist('lognormal','MOM',[m,s]);
% y = abs(icdf(pd,0.05)-u);
% end

% close all
% R_c = beta0;
% fun = @(x) (R_c - logninv(0.05,x,sqrt(log(0.1414^2+1))))^2;
% R_mu = fminsearch(fun,R_c);
% R_sig = sqrt(log(0.1414^2+1));
% pd = ERADist('lognormal','PAR',[R_mu,R_sig]);
% icdf(pd,0.05)
% 
% pi_pdf = [  ERADist('lognormal','PAR',[R_mu,R_sig]);
%             ERADist('GEV','PAR',[gev_para(1),gev_para(2),gev_para(3)])];


% % close all
% R_mu = 1.2262 * beta0;
% R_std = 0.1414*R_mu;
% pd = ERADist('lognormal','MOM',[R_mu,R_std]);
% icdf(pd,0.05)
% 
% pd = ERADist('lognormal','MOM',[6.46e6,8.79e5]);
% icdf(pd,0.05)
%
% fun = @(m) findlognmean(m,0.1414,beta0);
% R_mu = fminsearch(fun,beta0);
% R_std = 0.1414*R_mu;
% pd = ERADist('lognormal','MOM',[R_mu,R_std]);
% icdf(pd,0.05)