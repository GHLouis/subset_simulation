function [x_vehicle,D,D_max]=cal_D_max(K,V,W)
%displacement
%return a length-t vector of value of interest at x_section
L=1000;           %bridge length (m)
x_section=L/2;    %location of measurement on the bridge
t=3600;
dt = 10;
D=zeros(length(V),t);
x_vehicle=zeros(length(V),t);
EI = 1e10;
F = W*9.807;
for i=1:length(V)
    for j=((i-1)*dt+1):t
        % update each potential car's location
        x_vehicle(i,j)=x_vehicle(i,j)+V(i)*(j-(i-1)*dt);
        if and(x_vehicle(i,j)>=0,x_vehicle(i,j)<L)
            % update contribution based on potential car's new location
            D(i,j)=F(i)*K(i)*((x_section*L/(6*EI)*(L-x_vehicle(i,j))-x_section*(L-x_vehicle(i,j))^3/(6*L*EI)-x_section^3*(L-x_vehicle(i,j))/(6*L*EI))*(x_section<=x_vehicle(i,j))...
                +(-x_section^3*(L-x_vehicle(i,j))/(6*L*EI)+x_section/(6*EI*L)*((L-x_vehicle(i,j))*(L^2-(L-x_vehicle(i,j)).^2))+1/(6*EI)*(x_section-x_vehicle(i,j))^3)*(x_section>x_vehicle(i,j)));
        else
            % ignore this potential car for the rest of the time
            break
        end
    end
end
D_max=sum(D);
% figure(2)
% plot(D_max);
end
