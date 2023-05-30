Pf_SuS      = zeros([Np 1]);    % failure probability
delta_SuS   = zeros([Np 1]);    % cov of each run
b           = cell([Np 1]);     % intermediate failure levels
Pf          = cell([Np 1]);     % intermediate failure probabilities
b_sus       = cell([Np 1]);     % limit state function values
pf_sus      = cell([Np 1]);     % failure probabilities corresponding to b_sus

fprintf('Subset simulation running...');
tic
for i=1:Np
    [Pf_SuS(i), delta_SuS(i), b{i}, Pf{i}, b_sus{i}, pf_sus{i}, samplesU, ~] = SuS(N,p0,g,pi_pdf);
    clc
    progress = round(i/Np*100);
    fprintf('Subset simulation running..., Progress: %d%%\n', progress);
    elapsedTime = toc;
    disp(['Current Time: ' datestr(now) ', Elapsed Time: ' num2str(elapsedTime)])
end
time = toc;
% save the result
t = datestr(datetime('now'),'yyyymmdd');
save([probname,'_',t,'.mat'],'Pf_SuS', 'delta_SuS', 'b', 'Pf', 'b_sus',...
    'pf_sus','time','beta','Np','pf_ex','N','p0','figpath','probname','samplesU')