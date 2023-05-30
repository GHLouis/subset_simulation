function y = findlognmean(m,c,gev_para)
s = c*m;
pi_pdf = [  ERADist('lognormal','MOM',[m,s]);
    ERADist('GEV','PAR',[gev_para(1),gev_para(2),gev_para(3)])];
% correlation matrix
Rx = eye(2);   % independent case
% object with distribution information
pi_pdf = ERANataf(pi_pdf,Rx);    % if you want to include dependence
% limit state function and its gradient in the original space
g    = @(x) x(1)-x(2);
% Solve the optimization problem of the First Order Reliability Method
% FORM using MATLAB fmincon
[~, ~, ~,pf_fmc] = FORM_fmincon(g, pi_pdf);
y = abs(pf_fmc-1e-6);
end