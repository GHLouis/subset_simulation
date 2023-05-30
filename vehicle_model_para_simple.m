% 这个脚本用于绘制车辆荷载随机模型的参数的PDF
clc;clear;close all
path='figures\';
%% influence line function
% for practicle brigde, we use FEM to obtain the influence line function.
% every time we move the unit load and calculate the corresponse.

% # discretelize influence line function
L   = 35;                       % span
effect_type = 'moment';         % type of vehicle load effect
% effect_type = 'shearforce';
% effect_type = 'supportmoment';
switch effect_type
    case 'moment'
        % moment
        x   = [0;L/2;L];                % critical position
        y   = [0;1/4*L;0];              % critical ordinate of influence line function
        IL  = fit(x,y,'linearinterp');  % linear interpolation
    case 'shearforce'
        % shear force
        x   = [0;L];                % critical position
        y   = [1;0];              % critical ordinate of influence line function
        IL  = fit(x,y,'linearinterp');  % linear interpolation
    case 'supportmoment'
        L0 = L/2;
        fun = @(x) 1/(4*L0^2)*(L0-x).*x.*(L0-x+L0);
        IL = @(x) fun(x-L/2).*(x>=L/2) + fun(L/2-x).*(x<L/2);
end
dx  = 0.1;                      % step length for discretelization

% plot the influence line function
x       = 0:dx:L;
Inline  = IL(x);
figure; plot(x, Inline)
grid on
xlabel('Position of unit load (m)')
ylabel('Response (m)')
saveas(gcf,[path,'influence line function.png'])
savefig(gcf,[path,'influence line function'],'compact')

% test
xnew    = 0.3*L;
Inline(floor(xnew/dx)+1)
IL(xnew)

% % # maximum influence area under unit uniformly distributed loads
% A = zeros(L/dx,L/dx);
% Lc = dx:dx:L;
% Lc = Lc(:);
% if L <=100
%     for i = 1: L/dx                             % convoy length
%         lc = Lc(i);
%         for j = i: (L-lc)/dx+i                  % position of the first vehicle
%             x       = (0:dx:lc) + (j-i)*dx;
%             x       = x(:);
%             y       = Inline(floor(x/dx)+1);
%             A(i, j) = trapz(x, y);
%         end
%     end
%     mifA = [0; max(A,[],2)];
%     fun = @(x) x.*(2*L-x)/8;
%     figure
%     hold on
%     plot([0;Lc],mifA)
%     plot(Lc,fun(Lc))
%     grid on
%     xlabel('convoy length (m)')
%     ylabel('maximum influence area')
%     legend({'numerical solution','exact solution'},'Location','southeast')
%     saveas(gcf,[path,'maximum influence area.png'])
%     savefig(gcf,[path,'maximum influence area'],'compact')
% else
mifA = nan;
% end

% store result of influence line function
info.L              = L;
info.dx             = dx;
info.Inline         = Inline;
info.mifA           = mifA;
info.mInline        = max(Inline); % max of influence line function
%% vehicle: wheelbase and axle weight ratio
% accoring to zhang's book, table 4-6
wheelbase=[ 3       0       0       0       0;
    5       0       0       0       0;
    5       1.3     0       0       0;
    2.5     6       1.3     0       0;
    3.4     7.4     1.3     1.3     0;
    3.2     1.5     7       1.3     1.3];

wratio=[0.45	0.55	0       0       0       0;
    0.28	0.72	0       0       0       0;
    0.15	0.44	0.41	0       0       0;
    0.10	0.19	0.36	0.35	0       0;
    0.060	0.26	0.24	0.22	0.22	0;
    0.040	0.19	0.16	0.21	0.19	0.21];

% position of the resultant for relative to the front axle.
X   = wratio(:,2:end)*(cumsum(wheelbase,2)');
xc  = diag(X);

% store result of wheelbase and axle weight ratio
info.wheelbase      = wheelbase;
info.wratio         = wratio;
info.xc             = xc;
%% Gross vehicle weight
% see: "车辆荷载概率分布类型及其参数” (李 et al., 1997, p. 104)
% unit: ton

% pdw = makedist('Lognormal','mu',1.666697,'sigma',0.816272);
pdw         = cell(6,1);

pdw{1}      = makedist('norm','mu',150,'Sigma',50);
pdw{2}      = makedist('norm','mu',250,'Sigma',100);
pdw{3}      = makedist('norm','mu',350,'Sigma',200);
pdw{4}      = makedist('norm','mu',450,'Sigma',200);
pdw{5}      = makedist('norm','mu',550,'Sigma',200);
pdw{6}      = makedist('norm','mu',650,'Sigma',200);

info.pdw_new        = pdw;
%% inter-vehicle distance
% free headway, larger than 3 sec
pd_hf = makedist('Gamma',"a",0.904286,"b",1/0.039451);
% congestion headway, less than 3 sec
pd_hc = makedist('Gamma',"a",12.907330,"b",1/7.235810);

% 前车后轴与后车前轴之间的距离
% inter-vehicle distance = headway * speed;

pd_df = makedist('Lognormal','mu',4.827692,'sigma',1.115751); % 自由车距
pd_dc = makedist('Lognormal','mu',1.561165,'sigma',0.279707); % 拥堵车距

info.pd_hf = pd_hf;
info.pd_hc = pd_hc;
info.pd_df = pd_df;
info.pd_dc = pd_dc;
%% bridge design code
% I level load
qk = 10.5;                  % unit: kN/m
if info.L <= 5
    Pk = 270;               % unit: kN
elseif info.L > 5 && info.L < 50
    Pk = 2*(info.L + 130);  % unit: kN
else
    Pk = 360;               % unit: kN
end
M1_G    = info.L^2/8*qk;    % global loading
M1_L    = Pk*info.mInline;  % local loading
M1      = M1_G + M1_L;

% II level load
M2_G    = 0.75*M1_G;         % global loading
M2_L    = 0.75*M1_L;         % local loading
M2      = M2_G + M2_L;

info.M1  = M1;
info.M2  = M2;

save(['vehicle-load-model_',num2str(L),'m_',effect_type,'.mat'],'info');