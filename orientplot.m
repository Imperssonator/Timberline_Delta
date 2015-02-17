function out = orientplot(ORIENT,MAXCONF)

figure
% cutoff = 0.95;
% pcolor((MAXCONF>cutoff).*ORIENT + (MAXCONF<=cutoff).*(MAXCONF~=0).*-120 + (MAXCONF==0).*-180);

pcolor(ORIENT + (MAXCONF==0).*-180);
hc=colorbar; shading flat; axis equal; set(gca,'YDir','reverse');

out = 0;

end