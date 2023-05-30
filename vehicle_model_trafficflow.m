% define parameters of the random variables
d = 1000;                   % AADT, q=1000, 5000, 10000, 20000
% ratio of traffic volume per day (24 hour)
P_q = zeros(24,1) + 1/24;   % take the average, for example
q = mnrnd(d, P_q);
qcs = cumsum(q);

% ratio of vehicle type within each hour per day
P_cls = zeros(24,6) + 1/6;    % take the average, for example
% truck is not allowed to pass away in 7:9 and 17:20 in Guangzhou.
P_cls([7:9 17:20],2:end)    = 0;
P_cls([7:9 17:20],1)        = 1;

pd_cls      = cell(24,1);
for i=1:24
    pd_cls{i} = makedist('Multinomial','Probabilities',P_cls(i,:));
end

% hourly traffic volume
% to obtain the predefiend AADT, its randomness is not considered
pd_q   = makedist('Multinomial','Probabilities',P_q);

save('trafficflow.mat','d','pd_cls','pd_q','q','qcs')