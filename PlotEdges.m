function out = PlotEdges(XYZ,Edge)

figure
hold on

for i = 1:length(Edge)
    disp(i)
    p1 = XYZ(Edge(i,1),1:2);
    p2 = XYZ(Edge(i,2),1:2);
    plot([p1(1); p2(1)],[p1(2); p2(2)],'-b')
end

end