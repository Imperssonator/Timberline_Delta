function out = PickAmorphE()

Avg = -0.3; % Actually -4.9 but we're setting P3HT crystal HOMO as 0 (-4.6)
% In eV

SD = 0.05; % experimental error was 0.1 = 2sigma

out = randn()*SD+Avg;

end