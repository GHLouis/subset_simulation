function [x_vehicle,R,R_max]=cal_R_max(K,V,W)
%left end support reaction
%return a length-t vector of value of interest at x_section
L=1000;           %bridge length (m)
%x_section=0;    %location of measurement on the bridge
t=3600;
dt = 10;
R=zeros(length(V),t);
x_vehicle=zeros(length(V),t);
F = W*9.807;
for i=1:length(V)
    for j=((i-1)*dt+1):t
        % update each potential car's location
        x_vehicle(i,j)=x_vehicle(i,j)+V(i)*(j-(i-1)*dt);
        if and(x_vehicle(i,j)>=0,x_vehicle(i,j)<L)
            % update contribution based on potential car's new location
            R(i,j)=F(i)*K(i)*((L-x_vehicle(i,j))/L);
        else
            % ignore this potential car for the rest of the time
            break
        end
    end
end
R_max=sum(R);
% figure(3)
% plot(R_max);
end

% figure(1)
% plot(M_max,'r')
% hold on
% plot(R_max.*max(M_max)/max(R_max),'b--')
%  hold on
% plot(D_max.*max(M_max)/max(D_max),'y-.')
% figure(2)
% for i=1:360
%     hold on
%     plot(x_vehicle(i,:))
% end
% figure(3)
% plot(v)