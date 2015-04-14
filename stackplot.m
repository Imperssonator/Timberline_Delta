function out = stackplot(Stack)

%% Stackplot
% This function takes a pi stack matrix [start end] Nx2 and plots it as evenly
% spaced line segments


Colorset = [1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1; 0 0 0; 1 1 0; 0.8 0.8 0.8; 0.5 0 0; 0.6 0 0.6; 0.9 0 0.9];

for i = 1:size(Stack,3)
    for j = 1:size(Stack,1)
        if i<=length(Colorset)
            plot3([j j],[i i],[Stack(j,1,i),Stack(j,2,i)],'Color',Colorset(i,:))
        else
            plot3([j j],[i i],[Stack(j,1,i),Stack(j,2,i)],'Color',Colorset(mod(i,length(Colorset)),:))
        end
%     plot3([i,i], [j j], [Stack(i,1),Stack(i,2)],'Color',Colorset(j,:))
        hold on
    end
    drawnow
end

% rotate(gca,[0 0 1],45);
az = -11;
el = 50;
view(az, el);
end
