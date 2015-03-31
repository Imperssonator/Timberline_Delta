function out = stackplot(Stack,iter)

%% Stackplot
% This function takes a pi stack matrix [start end] Nx2 and plots it as evenly
% spaced line segments

Colorset = [1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0; 1 1 0; 0.8 0.8 0.8; 0.5 0 0; 0.6 0 0.6];
j = iter/10000;

for i = 1:size(Stack,1)
    plot3([i,i],[j,j],[Stack(i,1),Stack(i,2)],'Color',Colorset(iter/10000,:))
%     hold on
end

% rotate(gca,[0 0 1],45);
az = -11;
el = 50;
view(az, el);
end
