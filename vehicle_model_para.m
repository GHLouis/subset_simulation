% 这个脚本用于绘制车辆荷载随机模型的参数的PDF

% 车重W：只考虑总车重，用双峰高斯分布表示
% 车长：不考虑
% 车流量Q：daily traffic volume，用Weibull分布表示；同时，每个小时的交通量是变化的，用小时的平均值表示。
% 时间间距：小于4s，与车流量相关；大于4秒，与车流量无关，用移位指数分布
% 车速V：高斯分布
clear;close all
path='figures\';
%% 车重w
mu = [15;35];                    % Means
sigma = cat(3,[10],[20]); % Covariances
p = [0.6 0.4];                     % Mixing proportions
pdw = gmdistribution(mu,sigma,p);

x=random(pdw,1e3);
% dist1 = optimal_commom_dist(x);
% dist2 = optimal_kernel_dist(x);
dist3 = gmfit_best(x);

figure('Units','centimeters','Position',[18,7,10,8]);
x = linspace(0,60);
y = pdf(pdw,x');
plot(x,y,'k','linewidth',1.5);
xlabel('GVW(ton)');
ylabel('PDF');
set(gca,'Fontsize',10,'FontName','Times New Roman')
f_name='weight';
% exportgraphics(gcf,[path,f_name,'.pdf'])
% exportgraphics(gcf,[path,f_name,'.png'],'resolution',600)
% savefig(gcf,[path,f_name],'compact')

% 线性插值拟合车重分布的icdf函数
x = transpose(linspace(0,70,1e2));
y = cdf(pdw,x);
cdf_gmm=fit(x,y,'linearinterp');
icdf_gmm=fit(y,x,'linearinterp');
y1=cdf_gmm(x);
figure('Units','centimeters','Position',[18,7,14,10]);
subplot(2,2,1)
plot(x,y,'k.','linewidth',1.5);
hold on
plot(x,y1)
xlabel('GVW(ton)');
ylabel('CDF');
grid on
set(gca,'Fontsize',10,'FontName','Times New Roman')

subplot(2,2,2)
plot(x,-log(-log(y)),'k.','linewidth',1.5);
hold on
plot(x,-log(-log(y1)))
xlabel('GVW(ton)');
ylabel('$-ln(-ln(CDF))$','Interpreter','latex');
grid on
set(gca,'Fontsize',10,'FontName','Times New Roman')

subplot(2,2,3)
plot(y,x,'k.','linewidth',1.5);
hold on
plot(y1,x)
ylabel('GVW(ton)');
xlabel('CDF');
grid on
set(gca,'Fontsize',10,'FontName','Times New Roman')

subplot(2,2,4)
plot(-log(-log(y)),x,'k.','linewidth',1.5);
hold on
plot(-log(-log(y1)),x)
ylabel('GVW(ton)');
xlabel('$-ln(-ln(CDF))$','Interpreter','latex');
grid on
set(gca,'Fontsize',10,'FontName','Times New Roman')

f_name='weight_fitting';
% exportgraphics(gcf,[path,f_name,'.pdf'])
% exportgraphics(gcf,[path,f_name,'.png'],'resolution',600)
% savefig(gcf,[path,f_name],'compact')
%% 车流量q
% 当形状系数\xi<0的时候，GEV分布有上限：\mu-\sigma/\xi
% 现在取：\xi=-0.1，\mu=4800，\sigma=384，使得上限为8640
pdq = makedist('gev',-.1,384,4800);
x=random(pdq,1e3,1);
dist1 = optimal_commom_dist(x);
% dist2 = optimal_kernel_dist(x);
% dist3 = gmfit_best(x);

figure('Units','centimeters','Position',[18,7,10,8]);
x = linspace(3000,9000);
y = pdf(pdq,x);
plot(x,y,'k','linewidth',1.5);
xlabel('Traffic volume(vehs/day)');
ylabel('PDF');
set(gca,'Fontsize',10,'FontName','Times New Roman')
f_name='Traffic_volume';
% exportgraphics(gcf,[path,f_name,'.pdf'])
% exportgraphics(gcf,[path,f_name,'.png'],'resolution',600)
% savefig(gcf,[path,f_name],'compact')

% 车流量比例
Qratio = [zeros(1,6)+125 zeros(1,3)+175 zeros(1,7)+100 zeros(1,8)+150];
Qratio = Qratio./sum(Qratio);
figure('Units','centimeters','Position',[18,7,10,8]);
plot(Qratio,'k','linewidth',1.5);
xlabel('Time(h)');
ylabel('Ratio');
set(gca,'Fontsize',10,'FontName','Times New Roman')
f_name='Traffic_volume_ratio';
% exportgraphics(gcf,[path,f_name,'.pdf'])
% exportgraphics(gcf,[path,f_name,'.png'],'resolution',600)
% savefig(gcf,[path,f_name],'compact')
%% 车速
% 车速和车流量，用Green Shields模型。
Vm=20;
Qm=360;
% 假设车速和车流密度服从线性关系：即Green Shields模型。那么平均车速和车流量的关系为：
% Q=-4*Qm/Vm^2*V(V-Vm);
a=-4*Qm/Vm^2;
v2q = @(v) a.*v.*(v-Vm);
q2v_f = @(q) (a*Vm-sqrt((a*Vm)^2+4*a*q))/(2*a);
q2v_c = @(q) (a*Vm+sqrt((a*Vm)^2+4*a*q))/(2*a);

figure('Units','centimeters','Position',[18,7,10,8]);
v = linspace(0,Vm);
q = v2q(v);
plot(q,v,'k','linewidth',1.5);
hold on
plot([0  Qm],[Vm/2 Vm/2],'--','linewidth',1.5);
text(1/5*Qm,1/4*Vm,'Congenstion state','FontName','Times New Roman')
text(1/5*Qm,3/4*Vm,'Free flow state','FontName','Times New Roman')
ylabel('Average velocity(m/s)');
xlabel('Traffic volume(veh/hour)');
set(gca,'Fontsize',10,'FontName','Times New Roman')
f_name='velocity-volume';
% exportgraphics(gcf,[path,f_name,'.pdf'])
% exportgraphics(gcf,[path,f_name,'.png'],'resolution',600)
% savefig(gcf,[path,f_name],'compact')


% === 以下是草稿
% pd = makedist('norm',10,2);
% x=random(pd,1e3,1);
% dist1 = optimal_commom_dist(x);
% % dist2 = optimal_kernel_dist(x);
% % dist3 = gmfit_best(x);
% 
% figure('Position',[500,400,400,300],'Units','centimeters')
% x = linspace(0,24);
% y = pdf(pd,x);
% plot(x,y,'linewidth',1.5);
% xlabel('Velocity(m/s)');
% ylabel('PDF');
% set(gca,'Fontsize',10,'FontName','Times New Roman')
% f_name='Velocity';
% exportgraphics(gcf,[path,f_name,'.pdf'])
% exportgraphics(gcf,[path,f_name,'.png'],'resolution',600)
% savefig(gcf,[path,f_name],'compact')
%% headway
% 简单起见，认为所有车的headway都一样，是车流量的倒数。




% headway follows exponential distribution, and its mean is dependent on traffic volume.
% 一个重要的问题：如何区分大小headway？
% % 验证normalised headway
% close all
% Q=[55;155;255;355;455]./3600;
% mu=1./Q;
% xi=0:0.1:100;
% ti=0:0.01:5;
% figure
% f1=subplot(1,1,1);
% hold on
% figure
% f2=subplot(1,1,1);
% hold on
% for i=1:5
%     y=exppdf(xi,mu(i));
%     z=exppdf(ti,1);
%     plot(f1,xi,y,'linewidth',1.5);
%     plot(f2,ti,z,'linewidth',1.5);
% end
% legend(f1,{'55';'155';'255';'355';'455'})
% legend(f2,{'55';'155';'255';'355';'455'})
% 
% 
% Q=[20;40;60;80;120];
% 300./Q

% L=100;
% x = [0;L/2;L];
% y = [0;1/4*L;0];
% % f0=fit(x,y,'linearinterp');
% % IL = @(x) f0(x).*(x>=0&x<=L); %影响线函数
% IL = fit(x,y,'linearinterp');
L = 35;
x = [0;L/2;L];
y = [0;1/4*L;0];
IL = fit(x,y,'linearinterp');
dx = 0.1;
x = 0:dx:L;
Inline = IL(x);

% L       = 100;
% x       = [0;L/2;L];
% y       = [0;1/4*L;0];
% IL      = fit(x,y,'linearinterp');
% dx      = 0.1;
% x       = 0:dx:L;
% Inline  = IL(x);
% plot(IL)

% xnew    = 0.1;
% Inline(floor(xnew/dx)+1)
% IL(xnew)

A = zeros(L/dx,L/dx);
Lc = dx:dx:L;
for i = 1: L/dx         % 确定车队的长度
    lc = Lc(i);
    for j = i: (L-lc)/dx+i     % 确定车队领头车的位置
        x = (0:dx:lc) + (j-i)*dx;
        y = Inline(floor(x/dx)+1);
        A(i,j) = trapz(x,y);
    end
end

mifA = [0; max(A,[],2)];
% figure
% hold on
% plot(Lc(:),Y)
% 
% fun = @(x) x.*(2*L-x)/8;
% plot(Lc(:),fun(Lc))
% xlabel('车队长度')
% ylabel('最大影响线面积')
% legend({'数值解','解析解'},'Location','best')
wheelbase=[3	0	0	0	0;
5	0	0	0	0;
5	1.30000000000000	0	0	0;
2.50000000000000	6	1.30000000000000	0	0;
3.40000000000000	7.40000000000000	1.30000000000000	1.30000000000000	0;
3.20000000000000	1.50000000000000	7	1.30000000000000	1.30000000000000];

wratio=[0.450000000000000	0.550000000000000	0	0	0	0;
0.280000000000000	0.720000000000000	0	0	0	0;
0.150000000000000	0.440000000000000	0.410000000000000	0	0	0;
0.100000000000000	0.190000000000000	0.360000000000000	0.350000000000000	0	0;
0.0600000000000000	0.260000000000000	0.240000000000000	0.220000000000000	0.220000000000000	0;
0.0400000000000000	0.190000000000000	0.160000000000000	0.210000000000000	0.190000000000000	0.210000000000000];

info.Qratio         = Qratio;
info.IL             = IL;
info.pdw            = pdw;
info.pdq            = pdq;
info.icdf_gmm       = icdf_gmm;
info.q2v_f          = q2v_f;
info.q2v_c          = q2v_c;

% 影响线
info.L              = L;
info.dx             = dx;
info.Inline         = Inline;
info.mifA           = mifA;
info.mInline        = max(Inline);

info.wheelbase      = wheelbase;
info.wratio         = wratio;

%% 改用新的车重、车间距、车辆类别的分布
% 车重的分布:单位是吨
% 参见："车辆荷载概率分布类型及其参数” (李 et al., 1997, p. 104)
pdw         = cell(6,1);
pdw{1}      = makedist('Lognormal','mu',1.666697,'sigma',0.816272);
pdw{2}      = makedist('Lognormal','mu',1.666697,'sigma',0.816272);
pdw{3}      = makedist('Lognormal','mu',1.666697,'sigma',0.816272);
pdw{4}      = makedist('Lognormal','mu',1.666697,'sigma',0.816272);
pdw{5}      = makedist('Lognormal','mu',1.666697,'sigma',0.816272);
pdw{6}      = makedist('Lognormal','mu',1.666697,'sigma',0.816272);

pdw         = cell(6,1);
pdw{1}      = makedist('unif','lower',100,'upper',200);
pdw{2}      = makedist('unif','lower',200,'upper',300);
pdw{3}      = makedist('unif','lower',300,'upper',400);
pdw{4}      = makedist('unif','lower',400,'upper',500);
pdw{5}      = makedist('unif','lower',500,'upper',600);
pdw{6}      = makedist('unif','lower',600,'upper',700);


% % 24小时车流量的比例
% P_q         = zeros(24,1) + 1/24;   % take the average

% 24小时，每小时各种车型占比
P_cls       = zeros(1,6) + 1/6;    % take the average

pd_cls      = makedist('Multinomial','Probabilities',P_cls);
pd_df       = makedist('Lognormal','mu',4.827692,'sigma',1.115751); % 自由车距
pd_dc       = makedist('Lognormal','mu',1.561165,'sigma',0.279707); % 拥堵车距

info.pdw_new        = pdw;
info.pd_cls         = pd_cls;
info.pd_df          = pd_df;
info.pd_dc          = pd_dc;

save 'vehicle-load-model.mat' 'info';