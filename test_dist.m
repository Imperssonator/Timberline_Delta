function out = test_dist(SD,DPavg)

sample = ceil(randn(1000,1).*SD+DPavg);
%disp(length(sample))

out = sample;

% table = tabulate(sample);                   % [DP, count, mol frac]
% out = table(table(:,2)~=0 & table(:,1)>0,:);
% 
% DPn = Mn(sample);
% DPw = Mw(sample);
% PDI = DPw/DPn;

end