function out = test_dist(SD,DPavg)

sample = ceil(randn(10000,1).*SD+DPavg);
%disp(length(sample))

table = tabulate(sample);                   % [DP, count, mol frac]
out = table(table(:,2)~=0,:);

DPn = Mn(sample);
DPw = Mw(sample);
PDI = DPw/DPn;

end