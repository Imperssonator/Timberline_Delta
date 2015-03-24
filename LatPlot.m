function out = LatPlot(XYZ,bins)

figure
hold on

plot(XYZ(XYZ(:,3)==bins,1),XYZ(XYZ(:,3)==bins,2),'.b')
plot(XYZ(XYZ(:,3)~=bins,1),XYZ(XYZ(:,3)~=bins,2),'.g')

axis equal

end