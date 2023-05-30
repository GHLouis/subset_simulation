clear; close all; clc;
%
% obrien_example_1
% % 
% Dai_example_4
% Dai_example_5
% Dai_example_6
%%
figpath = 'figures\';
load('obrien example 1_20230510.mat')
probname = 'obrien example 1';
% xlabelname = '$Y$ - Moment(kNm)';
xlabelname = '$Y$ - System response';
thre_opc    = 'b';
% define parameters of the random variables
d      = 3000;          % number of dimensions
pi_pdf = repmat(ERADist('normal','PAR',[40,5]), d, 1);

% limit state function
% beta = 1e10;
beta0 = pi_pdf(1).icdf(exp(log(1-1/1000)/(250*d)));
post_process_ss_cv

b_oe1_ex = b_ex;
b_oe1_m = mean(b_SuS);
b_oe1_std = std(b_SuS);
b_oe1_cov = std(b_SuS)/mean(b_SuS);
b_oe1_error = b_er_100;

beta_oe1_ex = beta_ex;
beta_oe1_m = mean(beta_SuS);
beta_oe1_std = std(beta_SuS);
beta_oe1_cov = std(beta_SuS)/mean(beta_SuS);
beta_oe1_error = beta_er_100;

Nt_oe1 = Na;

gev_para_oe1 = gev_para_mr;
%%
load('Dai example 4_20230510.mat')
probname = 'Dai example 4';
% xlabelname = '$Y$ - Moment(kNm)';
xlabelname = '$Y$ - System response';
thre_opc    = 'b';
p = [0.98 0.02];
k =  [0.1;-0.2];
sigma = [265;1400];
mu = [3500; 3600];

% define parameters of the random variables
d0     = 2000;
d      = d0 * 2;
fun = @(x) p(1)*gevcdf(x,k(1),sigma(1),mu(1)) ...
    + p(2)*gevcdf(x,k(2),sigma(2),mu(2));

y = exp(log(1-1/1000)/(250*d0));
f = @(x) abs(fun(x) - y);
% limit state function
% beta = 1e10;
beta0 = fminsearch(f,mean(mu));
post_process_ss_cv

b_de4_ex = b_ex;
b_de4_m = mean(b_SuS);
b_de4_std = std(b_SuS);
b_de4_cov = std(b_SuS)/mean(b_SuS);
b_de4_error = b_er_100;

beta_de4_ex = beta_ex;
beta_de4_m = mean(beta_SuS);
beta_de4_std = std(beta_SuS);
beta_de4_cov = std(beta_SuS)/mean(beta_SuS);
beta_de4_error = beta_er_100;

Nt_de4 = Na;
gev_para_de4 = gev_para_mr;
%%
load('Dai example 5_20230510.mat')
probname = 'Dai example 5';
% xlabelname = '$Y$ - Moment(kNm)';
xlabelname = '$Y$ - System response';
thre_opc    = 'b';

p = [0.8 0.2];
k =  [-0.02;0.07];
sigma = [190;230];
mu = [360; 180];

% define parameters of the random variables
d0     = 2000;
d      = d0 * 2;
fun = @(x) p(1)*gevcdf(x,k(1),sigma(1),mu(1)) ...
    + p(2)*gevcdf(x,k(2),sigma(2),mu(2));

y = exp(log(1-1/1000)/(250*d0));
f = @(x) abs(fun(x) - y);
% limit state function
% beta = 1e10;
beta0 = fminsearch(f,mean(mu));

post_process_ss_cv

b_de5_ex = b_ex;
b_de5_m = mean(b_SuS);
b_de5_std = std(b_SuS);
b_de5_cov = std(b_SuS)/mean(b_SuS);
b_de5_error = b_er_100;

beta_de5_ex = beta_ex;
beta_de5_m = mean(beta_SuS);
beta_de5_std = std(beta_SuS);
beta_de5_cov = std(beta_SuS)/mean(beta_SuS);
beta_de5_error = beta_er_100;

Nt_de5 = Na;
gev_para_de5 = gev_para_mr;
%%
load('Dai example 6_20230510.mat')
probname = 'Dai example 6';
% xlabelname = '$Y$ - Moment(kNm)';
xlabelname = '$Y$ - System response';
thre_opc    = 'b';

p = [0.5 0.3 0.2];
k =  [-0.05;nan;nan]; % gev,norm,norm
sigma = [80;120;105];
mu = [60; 215;190];

% define parameters of the random variables
d0     = 3500;
d      = d0 * 3;
fun = @(x) p(1)*gevcdf(x,k(1),sigma(1),mu(1)) ...
            + p(2)*normcdf(x,mu(2),sigma(2))...
            + p(3)*normcdf(x,mu(3),sigma(3));
%         
%         fun = @(x) p(1)*gevcdf(x,k(1),sigma(1),mu(1)) ...
%             + p(2)*normcdf(x,sigma(2),mu(2))...
%             + p(3)*normcdf(x,sigma(3),mu(3));

y = exp(log(1-1/1000)/(250*d0));
f = @(x) abs(fun(x) - y);
% limit state function
% beta = 1e10;
beta0 = fminsearch(f,mean(mu));

post_process_ss_cv

b_de6_ex = b_ex;
b_de6_m = mean(b_SuS);
b_de6_std = std(b_SuS);
b_de6_cov = std(b_SuS)/mean(b_SuS);
b_de6_error = b_er_100;

beta_de6_ex = beta_ex;
beta_de6_m = mean(beta_SuS);
beta_de6_std = std(beta_SuS);
beta_de6_cov = std(beta_SuS)/mean(beta_SuS);
beta_de6_error = beta_er_100;

Nt_de6 = Na;
gev_para_de6 = gev_para_mr;
%% 写入texfile
probname = 'Dai main';
fid = fopen(['../report\',probname,'.tex'], 'w');
% fprintf(fid, '\\documentclass{article}\n');
% fprintf(fid, '\\begin{document}\n');
fprintf(fid, '\\begin{table}[tbph]\n');
fprintf(fid, ['\\caption{Characteristic value of ',probname,' (%d runs)}\n'],Np);
fprintf(fid, '\\begin{centering}\n');
fprintf(fid, '\\begin{tabular}{cccccc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Example & Exact value & Mean & STD & COV & Relative bias (\\%%)\\tabularnewline\n');
fprintf(fid, '\\midrule\n');
fprintf(fid, 'Obrien example 1 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',b_oe1_ex, b_oe1_m,b_oe1_std,b_oe1_cov,b_oe1_error);
fprintf(fid, 'Dai example 4 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',b_de4_ex, b_de4_m,b_de4_std,b_de4_cov,b_de4_error);
fprintf(fid, 'Dai example 5 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',b_de5_ex, b_de5_m,b_de5_std,b_de5_cov,b_de5_error);
fprintf(fid, 'Dai example 6 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',b_de6_ex, b_de6_m,b_de6_std,b_de6_cov,b_de6_error);
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\par\\end{centering}\n');
fprintf(fid, '\\end{table}\n\n');



fprintf(fid, '\\begin{table}[tbph]\n');
fprintf(fid, ['\\caption{Estimated parameters of GEV distribution of ',probname,' (%d runs)}\n'],Np);
fprintf(fid, '\\begin{centering}\n');
fprintf(fid, '\\begin{tabular}{cccc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, '\\multirow{2}{*}{Example} & \\multicolumn{3}{c}{Parameters}\\tabularnewline\n');
fprintf(fid, ' & $\\xi$ & $\\sigma$ & $\\mu$\\tabularnewline\n');
fprintf(fid, '\\midrule\n');
fprintf(fid, 'Obrien example 1  & %.4f & %.2f & %.2f    \\tabularnewline\n',gev_para_oe1(1), gev_para_oe1(2),gev_para_oe1(3));
fprintf(fid, 'Dai example 4     & %.4f & %.2f & %.2f    \\tabularnewline\n',gev_para_de4(1), gev_para_de4(2),gev_para_de4(3));
fprintf(fid, 'Dai example 5     & %.4f & %.2f & %.2f    \\tabularnewline\n',gev_para_de5(1), gev_para_de5(2),gev_para_de5(3));
fprintf(fid, 'Dai example 6     & %.4f & %.2f & %.2f    \\tabularnewline\n',gev_para_de6(1), gev_para_de6(2),gev_para_de6(3));
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\par\\end{centering}\n');
fprintf(fid, '\\end{table}\n\n');


fprintf(fid, '\\begin{table}[tbph]\n');
fprintf(fid, ['\\caption{Reliability index of ',probname,' (%d runs)}\n'],Np);
fprintf(fid, '\\begin{centering}\n');
fprintf(fid, '\\begin{tabular}{cccccc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Example & Exact value & Mean & STD & COV & Relative bias (\\%%)\\tabularnewline\n');
fprintf(fid, '\\midrule\n');
fprintf(fid, 'Obrien example 1 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',beta_oe1_ex, beta_oe1_m,beta_oe1_std,beta_oe1_cov,beta_oe1_error);
fprintf(fid, 'Dai example 4 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',beta_de4_ex, beta_de4_m,beta_de4_std,beta_de4_cov,beta_de4_error);
fprintf(fid, 'Dai example 5 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',beta_de5_ex, beta_de5_m,beta_de5_std,beta_de5_cov,beta_de5_error);
fprintf(fid, 'Dai example 6 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',beta_de6_ex, beta_de6_m,beta_de6_std,beta_de6_cov,beta_de6_error);
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\par\\end{centering}\n');
fprintf(fid, '\\end{table}\n\n');

fclose(fid);