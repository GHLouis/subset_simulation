function [x_vehicle,M,M_max]=cal_M_max(K,V,W)

meth=1;

switch meth
    case 1
        % METHOD 1
        %moment
        %return a length-t vector of value of interest at x_section
        L=1000;           %bridge length (m)
        x_section=L/2;    %location of measurement on the bridge
        IL=@(x) 0.5*x.*(x>=0&x<=x_section)+(0.25*L-0.5*(x-x_section)).*(x>x_section&x<L);
        t=3600;
        dt = 10;
        M=zeros(length(V),t);
        x_vehicle=zeros(length(V),t);
        F = W*9.807;
        for i=1:length(V)
            for j=((i-1)*dt+1):t
                % update each potential car's location
                x_vehicle(i,j)=x_vehicle(i,j)+V(i)*(j-(i-1)*dt);
                if and(x_vehicle(i,j)>=0,x_vehicle(i,j)<L)
                    % update moment contribution based on potential car's new location
                                        M(i,j)=F(i)*K(i)*(((L-x_section)*x_vehicle(i,j)/L)*(x_vehicle(i,j)<=x_section)+(x_section*(L-x_vehicle(i,j))/L)*(x_vehicle(i,j)>x_section));
%                     M(i,j)=F(i)*K(i)*IL(x_vehicle(i,j));
                else
                    % ignore this potential car for the rest of the time
                    break
                end
            end
        end
        M_max=max(sum(M));
        % figure(1)
        % plot(M_max);
        % figure(4)
        % plot(x_vehicle(2,:))

    case 2

        % METHOD 2
        % define the influence line
        L=1000;           %bridge length (m)
        x_section=L/2;    %location of measurement on the bridge
        IL=@(x) 0.5*x.*(x>=0&x<=x_section)+(0.25*L-0.5*(x-x_section)).*(x>x_section&x<L);
        M = zeros(360,3600);
        K_rep=repmat(K,[1 3600]);
        W_rep=repmat(W,[1 3600]);


        % figure
        % fplot(IL,[-100,1200])
        t=1:1:3600;
        % moving time of vehicle
        ch=0:10:3590;
        ch=ch(:);
        T=t-ch;
        % position of vehicle
        x_vehicle=T.*V;
        % find the index of effective x_vehicle
        id = x_vehicle>=0&x_vehicle< L;
        sum(id(:));
        % calculate the load effect
        M(id)=9.807*K_rep(id).*W_rep(id).*IL(x_vehicle(id));
        % obtain the max of each time step
        M_max=max(sum(M));
end
end
