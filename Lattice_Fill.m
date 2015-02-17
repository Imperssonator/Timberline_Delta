function [MonList, Edges] = Lattice_Fill(MS,pw,ph,lsx,lsy)

%% Lattice Fill
% MS is a binned microstructure where each pixel in a bin represents its
% angle off the x-axis

% pw = pixel width
% ph = pixel height
% lsx = lattice spacing x
% lsy = lattice spacing y

[h, w, bins] = size(MS);
h = h*ph;
w = w*pw;

%% Generate Monomer Grid

for i = 1:bins
    Img_Corners = [0 0;...
                   w 0;...
                   0 h;...
                   w h];
    MSi = MS(:,:,i);
    Bin_Angle = MSi(find(MSi,1))/180*pi;
    Img_Corners_Rot = Rotate_List(Img_Corners,Bin_Angle);
    
    [Top_Corner max_y] = max(Img_Corners_Rot(:,2));
    [Bottom_Corner min_y] = min(Img_Corners_Rot(:,2));
    [Right_Corner max_x] = max(Img_Corners_Rot(:,1));
    [Left_Corner min_x] = min(Img_Corners_Rot(:,1));
    Left_Side = min_x-1;        % Left side of a line segment that will completely span the diamond
    Right_Side = max_x+1;
    
    xbox = [Img_Corners_Rot(:,1)' Img_Corners_Rot(1,1)];
    ybox = [Img_Corners_Rot(:,2)' Img_Corners_Rot(1,2)];
    
    for y = min_y:lsy:max_y
        
        % Some kind of code here that finds the intersection points
        
    end
    
    
end

end

function Rotated = Rotate_List(PtList,angle)

ROT = [cos(angle) -sin(angle);...
       sin(angle) cos(angle)];
Rotated = [];
   for i = 1:size(PtList,1)
       Rotated(i,:) = ROT*PtList(i,:)';
   end
end