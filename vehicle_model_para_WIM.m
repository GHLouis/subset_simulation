% 这个脚本用于绘制车辆荷载随机模型的参数的PDF
clc;clear;close all
path='figures\';
p_name = 'G107';
%% influence line function
% for practicle brigde, we use FEM to obtain the influence line function.
% every time we move the unit load and calculate the corresponse.

% # discretelize influence line function
LL0   = [15;35];                       % span
% type of vehicle load effect
effect_type0 = {'moment';'shearforce';'supportmoment'};

for NL = 1:2
    for NE = 1:3
        L = LL0(NL);
        effect_type = effect_type0{NE};
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
        % saveas(gcf,[path,'influence line function.png'])
        % savefig(gcf,[path,'influence line function'],'compact')

        % test
        xnew    = 0.3*L;
        Inline(floor(xnew/dx)+1)
        IL(xnew)


        mifA = nan;

        % store result of influence line function
        info.L              = L;
        info.dx             = dx;
        info.Inline         = Inline;
        info.mifA           = mifA;
        info.mInline        = max(Inline); % max of influence line function
        %% vehicle: wheelbase and axle weight ratio
        wheel = [ 2.6       0       0       0       0;
            4.5       0       0       0       0;
            4.0       0       0       0       0;
            3.5       1.3     0       0       0;
            2     5.5       1.3     0       0;
            3     6     1.3     1.3     0;
            3     1.3     6.5       1.3     1.3];

        wr =[0.44	0.56	0       0       0       0;
            0.46	0.54	0       0       0       0;
            0.35	0.65	0       0       0       0;
            0.20	0.40	0.40	0       0       0;
            0.15	0.31	0.27	0.27	0       0;
            0.10	0.30	0.20	0.20	0.20	0;
            0.06	0.20	0.20	0.18	0.18	0.18];
        % Gross vehicle weight

        load([p_name,'_GVW fitting.mat'])
        for j = 1:2 % Direction
            for k = (j-1)*3 + (1:3) % Lane
                probname = [p_name,'_trafficflow_','Dr_',num2str(j),'_Lane_',num2str(k)];
                load([probname,'.mat'])

                pdw = []; % unit: kN
                CLS = [];
                wheelbase = [];
                wratio = [];
                Numaxle = [];
                for m = 1:7
                    NumComp = gm_model{m, 1}.NumComponents;
                    p = gm_model{m, 1}.ComponentProportion;
                    mu = gm_model{m, 1}.mu;
                    sigma = sqrt(gm_model{m, 1}.Sigma(:));
                    for kk = 1:NumComp
                        pd = makedist('norm','mu',mu(kk),'Sigma',sigma(kk));
                        pdw = [pdw; pd];
                    end
                    CLS = [CLS cls(:,m)*p];
                    wheelbase = [wheelbase; repmat(wheel(m,:),[NumComp 1])];
                    wratio = [wratio; repmat(wr(m,:),[NumComp 1])];
                    if m <=3
                        Numaxle = [Numaxle;zeros(NumComp,1) + 2];
                    else
                        Numaxle = [Numaxle;zeros(NumComp,1) + m-1];
                    end
                end

                pd_cls  = cell(24,1);
                for i=1:24
                    pd_cls{i} = makedist('Multinomial','Probabilities',CLS(i,:));
                end

                info.pdw_new  = pdw;
                info.Npdw     = length(pdw);
                info.d        = AADT;
                info.q        = round(AADT*q);
                info.qcs      = cumsum(info.q);
                info.one      = ones(sum(info.q),1);
                info.pd_cls   = pd_cls;
                info.Numaxle   = Numaxle;
                %% inter-vehicle distance
                % info.pd_d = pdd;
                info.pd_d = pdvgd;
                %% vehicle: wheelbase and axle weight ratio
                % wheelbase = [ 2.6       0       0       0       0;
                %             4.5       0       0       0       0;
                %             4.0       0       0       0       0;
                %             3.5       1.3     0       0       0;
                %             2     5.5       1.3     0       0;
                %             3     6     1.3     1.3     0;
                %             3     1.3     6.5       1.3     1.3];
                %
                % wratio =[0.44	0.56	0       0       0       0;
                %         0.46	0.54	0       0       0       0;
                %         0.35	0.65	0       0       0       0;
                %         0.20	0.40	0.40	0       0       0;
                %         0.15	0.31	0.27	0.27	0       0;
                %         0.10	0.30	0.20	0.20	0.20	0;
                %         0.06	0.20	0.20	0.18	0.18	0.18];

                % position of the resultant for relative to the front axle.
                X   = wratio(:,2:end)*(cumsum(wheelbase,2)');
                xc  = diag(X);

                % store result of wheelbase and axle weight ratio
                info.wheelbase      = wheelbase;
                info.wratio         = wratio;
                info.xc             = xc;
                info.VL             = sum(wheelbase,2);
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

                save([p_name,'_trafficflow_','Dr_',num2str(j),'_Lane_',num2str(k),'_WIM-vehicle-load-model_',num2str(L),'m_',effect_type,'.mat'],'info');
            end
        end
    end
end