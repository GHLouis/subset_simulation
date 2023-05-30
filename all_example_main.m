Dai_main
Practicle_24_main
Practicle_main
Practicle_WIM_main

%曲线标签
repalce_legend

%%
folder_path = pwd;
folder = [folder_path,'\figures'];  % 替换为实际文件夹路径
cd(folder)
% 手动执行 compile_all_tex.bat 文件
system('compile_all_tex.bat')
cd(folder_path)
close all
clc