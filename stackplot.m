function out = stackplot(Stack)

%% Stackplot
% This function takes a pi stack matrix [start end] Nx2 and plots it as evenly
% spaced line segments

for i = 1:size(Stack,1)
    plot([i,i],[Stack(i,1),Stack(i,2)],'-b')
%     hold on
end

end
