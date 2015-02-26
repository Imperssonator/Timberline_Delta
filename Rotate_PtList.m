function Rotated = Rotate_PtList(PtList,angle)

ROT = [cos(angle) -sin(angle);...
       sin(angle) cos(angle)];
Rotated = zeros(size(PtList,1),2);
   for i = 1:size(PtList,1)
       Rotated(i,:) = ROT*PtList(i,:)';
   end
end