function [M,idM] = maxLE_axle_gd(vehicle,info)
% modified on 20230519
% INPUT====================================================================
% vehicle:  vehicle type, inter-vehicle distance, weight, and etc. matrix
% info:     vehicel load model. strcture body
% OUTPUT===================================================================
% M:        target load effect. scalar
% =========================================================================
n           =   size(vehicle,1); % num of vehicles
% dx          =   info.dx;       % step for calculation
dx          =   0.1; % step for calculation

% initial position of front axle of each vechile
% assume vehicle length is 8m, center of weight is at the middle of
% vehicle.
x0 = [0;-cumsum(vehicle(1:end-1,4))];
x0(2:end) = x0(2:end) - cumsum(vehicle(1:end-1,2));
n_dx = ceil((-x0(end)+info.L)/dx); % No. of vehicle moving steps
% M = zeros(n_dx,n);
% m = zeros(n_dx,1);
% 
% % obtain load effect vehicle by vehicle
% for i = 1 : n
%     if vehicle(i,3) > 300 % 超过30吨的车，才有必要算
%         j = -ceil(x0(i)) + (dx:info.L);
%         x_vehicle = x0(i) + j*dx;
% 
%         %         % vehicle by vehicle
%         %         m(j) = vehicle(i,3)*info.Inline(floor(x_vehicle/info.dx)+1);
% 
%         % current vehicle, axle by axle
%         x_k = x_vehicle(:) - cumsum(vehicle(i,11:11+vehicle(i,1)-1));
%         x_k = min(max(x_k,0),info.L);
%         m(j) = info.Inline(floor(x_k/info.dx)+1) * transpose(vehicle(i,5:5+vehicle(i,1)-1));
% 
%         M(:,i) = m;
%         m(j)=0;
%     end
% end
% [M,idM] = max(sum(M,2)); 
M = zeros(n_dx,1);
m = zeros(n_dx,1);

% obtain load effect vehicle by vehicle
for i = 1 : n
    if vehicle(i,3) > 300 % 超过30吨的车，才有必要算
        j = (-ceil(x0(i)) + (dx:dx:info.L))/dx;
        j = round(j);
        x_vehicle = x0(i) + j*dx;

        %         % vehicle by vehicle
        %         m(j) = vehicle(i,3)*info.Inline(floor(x_vehicle/info.dx)+1);

        % current vehicle, axle by axle
        x_k = x_vehicle(:) - cumsum(vehicle(i,11:11+vehicle(i,1)-1));
        x_k = min(max(x_k,0),info.L);
        m(j) = info.Inline(floor(x_k/info.dx)+1) * transpose(vehicle(i,5:5+vehicle(i,1)-1));

        M(j) = M(j) + m(j);
        m(j)=0;
    end
end
[M,idM] = max(M);
%% postprocess
x = x0 + idM*dx; % most unfavorable position of the queue
y = vehicle(:,3);
y0 = zeros(size(y));
ex = find(x>=0 & x<=info.L);
figure
p1 = subplot(2,1,1);
hold on
f = plot(x,y,'Marker','diamond','LineStyle','none');
ax = gca;

xlim([0 info.L])
ylim([0 ax.YLim(2)])

ax = gca;
axpos = ax.Position;

% renderer = get(gcf, 'Renderer');
n = length(ex);
% vehicle by vehicle 
for i = 1:n
    % axle by axle  
        x_vehicle = x(ex(i));
        x_k = x_vehicle - cumsum(vehicle(ex(i),11:11+vehicle(ex(i),1)-1));
        x_k = min(max(x_k,0),info.L);
        y_k = vehicle(ex(i),5:5+vehicle(ex(i),1)-1);
        for j = 1:vehicle(ex(i),1)
            arx = [x_k(j)  x_k(j)];
            ary = [y_k(j) 0];
            u = axpos(1) + (arx - ax.XLim(1)) / diff(ax.XLim) * axpos(3);
            v = axpos(2) + (ary - ax.YLim(1)) / diff(ax.YLim) * axpos(4);
            annotation('arrow',u,v);
        end
end
delete(f)
box('on')
xlabel('Position of vehicle load (m)','Interpreter','latex')
ylabel('Axle weight (kN)','Interpreter','latex')
set(gca,"TickLabelInterpreter","latex")

p2 = subplot(2,1,2);
hold on
x_line = 0:info.dx:info.L;
x_line = x_line(:);
plot(x_line,info.Inline)
xlim([0 info.L])

% vehicle by vehicle 
for i = 1:n
    % axle by axle  
        x_vehicle = x(ex(i));
        x_k = x_vehicle - cumsum(vehicle(ex(i),11:11+vehicle(ex(i),1)-1));
        x_k = min(max(x_k,0),info.L);
        for j = 1:vehicle(ex(i),1)
            arx = [x_k(j)  x_k(j)];
            ary = [info.Inline(floor(x_k(j)/info.dx)+1) 0];
            plot(arx,ary,'k')
        end
end
set(gca,"TickLabelInterpreter","latex")
xlabel('Position of unit load (m)','Interpreter','latex')
ylabel('Influnce function','Interpreter','latex')
box('on')
% M = y_k * info.Inline(floor(x_k/info.dx)+1);
end