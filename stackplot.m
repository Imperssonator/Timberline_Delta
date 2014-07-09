function out = stackplot(Stack)

%% Stackplot
% This function takes a pi stack matrix [start end] Nx2 and plots it as evenly
% spaced line segments
figure
hold on

for i = 1:length(Stack)
    plot([i,i],[Stack(i,1),Stack(i,2)],'-b')
end

end
