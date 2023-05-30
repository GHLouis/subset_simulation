clear; close all; clc;
%
% % pf_ex = 1 - exp(log(1-1/100)/(365*24))
% Practicle_example_24_1
% Practicle_example_24_2
% Practicle_example_24_3
% Practicle_example_24_4
% Practicle_example_24_5
%%
figpath = 'figures\';
load('Practicle example 24 1_20230503.mat')
probname = 'Practicle example 24 1';
xlabelname = 'Moment (kNm)';
% xlabelname = 'Shear force (kNm)';
thre_opc    = 'a';
post_process_ss_cv

b_pe1_ex = b_ex;
b_pe1_m = mean(b_SuS);
b_pe1_std = std(b_SuS);
b_pe1_cov = std(b_SuS)/mean(b_SuS);
b_pe1_error = b_er_100;

beta_pe1_ex = beta_ex;
beta_pe1_m = mean(beta_SuS);
beta_pe1_std = std(beta_SuS);
beta_pe1_cov = std(beta_SuS)/mean(beta_SuS);
beta_pe1_error = beta_er_100;

Nt_pe1 = Na;
gev_para_pe1 = gev_para_mr;
%%
figpath = 'figures\';
load('Practicle example 24 2_20230503.mat')
probname = 'Practicle example 24 2';
% xlabelname = 'Moment (kNm)';
xlabelname = 'Shear force (kNm)';
thre_opc    = 'a';
post_process_ss_cv

b_pe2_ex = b_ex;
b_pe2_m = mean(b_SuS);
b_pe2_std = std(b_SuS);
b_pe2_cov = std(b_SuS)/mean(b_SuS);
b_pe2_error = b_er_100;

beta_pe2_ex = beta_ex;
beta_pe2_m = mean(beta_SuS);
beta_pe2_std = std(beta_SuS);
beta_pe2_cov = std(beta_SuS)/mean(beta_SuS);
beta_pe2_error = beta_er_100;

Nt_pe2 = Na;
gev_para_pe2 = gev_para_mr;
%%
figpath = 'figures\';
load('Practicle example 24 3_20230503.mat')
probname = 'Practicle example 24 3';
xlabelname = 'Moment (kNm)';
% xlabelname = 'Shear force (kNm)';
thre_opc    = 'a';
post_process_ss_cv

b_pe3_ex = b_ex;
b_pe3_m = mean(b_SuS);
b_pe3_std = std(b_SuS);
b_pe3_cov = std(b_SuS)/mean(b_SuS);
b_pe3_error = b_er_100;

beta_pe3_ex = beta_ex;
beta_pe3_m = mean(beta_SuS);
beta_pe3_std = std(beta_SuS);
beta_pe3_cov = std(beta_SuS)/mean(beta_SuS);
beta_pe3_error = beta_er_100;

Nt_pe3 = Na;
gev_para_pe3 = gev_para_mr;
%%
figpath = 'figures\';
load('Practicle example 24 4_20230503.mat')
probname = 'Practicle example 24 4';
% xlabelname = 'Moment (kNm)';
xlabelname = 'Shear force (kNm)';
thre_opc    = 'a';
post_process_ss_cv

b_pe4_ex = b_ex;
b_pe4_m = mean(b_SuS);
b_pe4_std = std(b_SuS);
b_pe4_cov = std(b_SuS)/mean(b_SuS);
b_pe4_error = b_er_100;

beta_pe4_ex = beta_ex;
beta_pe4_m = mean(beta_SuS);
beta_pe4_std = std(beta_SuS);
beta_pe4_cov = std(beta_SuS)/mean(beta_SuS);
beta_pe4_error = beta_er_100;

Nt_pe4 = Na;
gev_para_pe4 = gev_para_mr;
%%
figpath = 'figures\';
load('Practicle example 24 5_20230503.mat')
probname = 'Practicle example 24 5';
xlabelname = 'Moment (kNm)';
% xlabelname = 'Shear force (kNm)';
thre_opc    = 'a';
post_process_ss_cv

b_pe5_ex = b_ex;
b_pe5_m = mean(b_SuS);
b_pe5_std = std(b_SuS);
b_pe5_cov = std(b_SuS)/mean(b_SuS);
b_pe5_error = b_er_100;

beta_pe5_ex = beta_ex;
beta_pe5_m = mean(beta_SuS);
beta_pe5_std = std(beta_SuS);
beta_pe5_cov = std(beta_SuS)/mean(beta_SuS);
beta_pe5_error = beta_er_100;

Nt_pe5 = Na;
gev_para_pe5 = gev_para_mr;
%% 写入texfile
probname = 'Practicle main 24';
fid = fopen(['../report\',probname,'.tex'], 'w','n','UTF-8');
% fprintf(fid, '\\documentclass{article}\n');
% fprintf(fid, '\\begin{document}\n');
fprintf(fid, '\\begin{table}[tbph]\n');
fprintf(fid, ['\\caption{Characteristic value of ',probname,' (%d runs)}\n'],Np);
fprintf(fid, '\\begin{centering}\n');
fprintf(fid, '\\begin{tabular}{cccccc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Example & Exact value & Mean & STD & COV & Relative bias (\\%%)\\tabularnewline\n');
fprintf(fid, '\\midrule\n');
fprintf(fid, '1 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',b_pe1_ex, b_pe1_m,b_pe1_std,b_pe1_cov,b_pe1_error);
fprintf(fid, '2 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',b_pe2_ex, b_pe2_m,b_pe2_std,b_pe2_cov,b_pe2_error);
fprintf(fid, '3 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',b_pe3_ex, b_pe3_m,b_pe3_std,b_pe3_cov,b_pe3_error);
fprintf(fid, '4 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',b_pe4_ex, b_pe4_m,b_pe4_std,b_pe4_cov,b_pe4_error);
fprintf(fid, '5 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',b_pe5_ex, b_pe5_m,b_pe5_std,b_pe5_cov,b_pe5_error);
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
fprintf(fid, '1  & %.4f & %.2f & %.2f    \\tabularnewline\n',gev_para_pe1(1), gev_para_pe1(2),gev_para_pe1(3));
fprintf(fid, '2  & %.4f & %.2f & %.2f    \\tabularnewline\n',gev_para_pe2(1), gev_para_pe2(2),gev_para_pe2(3));
fprintf(fid, '3  & %.4f & %.2f & %.2f    \\tabularnewline\n',gev_para_pe3(1), gev_para_pe3(2),gev_para_pe3(3));
fprintf(fid, '4  & %.4f & %.2f & %.2f    \\tabularnewline\n',gev_para_pe4(1), gev_para_pe4(2),gev_para_pe4(3));
fprintf(fid, '5  & %.4f & %.2f & %.2f    \\tabularnewline\n',gev_para_pe5(1), gev_para_pe5(2),gev_para_pe5(3));
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\par\\end{centering}\n');
fprintf(fid, '\\end{table}\n\n');

fprintf(fid, '\\begin{table}[tbph]\n');
fprintf(fid, ['\\caption{Reliability inpex of ',probname,' (%d runs)}\n'],Np);
fprintf(fid, '\\begin{centering}\n');
fprintf(fid, '\\begin{tabular}{cccccc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Example & Exact value & Mean & STD & COV & Relative bias (\\%%)\\tabularnewline\n');
fprintf(fid, '\\midrule\n');
fprintf(fid, '1 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',beta_pe1_ex, beta_pe1_m,beta_pe1_std,beta_pe1_cov,beta_pe1_error);
fprintf(fid, '2 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',beta_pe2_ex, beta_pe2_m,beta_pe2_std,beta_pe2_cov,beta_pe2_error);
fprintf(fid, '3 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',beta_pe3_ex, beta_pe3_m,beta_pe3_std,beta_pe3_cov,beta_pe3_error);
fprintf(fid, '4 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',beta_pe4_ex, beta_pe4_m,beta_pe4_std,beta_pe4_cov,beta_pe4_error);
fprintf(fid, '5 & %.2f & %.2f & %.2f & %.2f & %.2f \\tabularnewline\n',beta_pe5_ex, beta_pe5_m,beta_pe5_std,beta_pe5_cov,beta_pe5_error);
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\par\\end{centering}\n');
fprintf(fid, '\\end{table}\n\n');
fclose(fid);