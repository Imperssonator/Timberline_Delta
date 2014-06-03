function [DPn,DPw,PDI] = test_dist(SD,DPavg)

sample = ceil(randn(10000,1).*SD+DPavg);
disp(length(sample))

DPn = Mn(sample);
DPw = Mw(sample);
PDI = DPw/DPn;

end