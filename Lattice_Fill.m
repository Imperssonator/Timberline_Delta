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
    
    [max_y ,Top_Corner] = max(Img_Corners_Rot(:,2));
    [min_y ,Bottom_Corner] = min(Img_Corners_Rot(:,2));
    [max_x ,Right_Corner] = max(Img_Corners_Rot(:,1));
    [min_x ,Left_Corner] = min(Img_Corners_Rot(:,1));
    Left_Side = min_x-1;        % Left side of a line segment that will completely span the diamond
    Right_Side = max_x+1;
    
    % Equations for the lines of the box
    X_bottomright = [Img_Corners_Rot(Bottom_Corner,1) Img_Corners_Rot(Right_Corner,1)];
    Y_bottomright = [Img_Corners_Rot(Bottom_Corner,2) Img_Corners_Rot(Right_Corner,2)];
    line_bottomright = polyfit(X_bottomright, Y_bottomright,1);
    
    X_topright = [Img_Corners_Rot(Top_Corner,1) Img_Corners_Rot(Right_Corner,1)];
    Y_topright = [Img_Corners_Rot(Top_Corner,2) Img_Corners_Rot(Right_Corner,2)];
    line_topright = polyfit(X_topright, Y_topright,1);
    
    X_topleft = [Img_Corners_Rot(Top_Corner,1) Img_Corners_Rot(Left_Corner,1)];
    Y_topleft = [Img_Corners_Rot(Top_Corner,2) Img_Corners_Rot(Left_Corner,2)];
    line_topleft = polyfit(X_topleft, Y_topleft,1);
    
    X_bottomleft = [Img_Corners_Rot(Bottom_Corner,1) Img_Corners_Rot(Left_Corner,1)];
    Y_bottomleft = [Img_Corners_Rot(Bottom_Corner,2) Img_Corners_Rot(Left_Corner,2)];
    line_bottomleft = polyfit(X_bottomleft, Y_bottomleft,1);
    
    xbox = [Img_Corners_Rot(:,1)' Img_Corners_Rot(1,1)];
    ybox = [Img_Corners_Rot(:,2)' Img_Corners_Rot(1,2)];
    
    Intersect = [Img_Corners_Rot(Bottom_Corner,:); Img_Corners_Rot(Top_Corner,:)];
    
    
    for y = min_y+lsy:lsy:max_y-lsy
        line_1 = [0 1 y]; % ax + by = c
        
        % Right intersection
        if y<=Img_Corners_Rot(Right_Corner,2)
            line_2 = [line_bottomright(1) -1 -1*line_bottomright(2)];
        else
            line_2 = [line_topright(1) -1 -1*line_topright(2)];
        end
        
        new_intersect1 = [line_1(1:2); line_2(1:2)]\[line_1(3);line_2(3)];
        
        % Left intersection
        if y <= Img_Corners_Rot(Left_Corner,2)
            line_2 = [line_bottomleft(1) -1 -1*line_bottomleft(2)];
        else
            line_2 = [line_topleft(1) -1 -1*line_topleft(2)];
        end
        
        new_intersect2 = [line_1(1:2); line_2(1:2)]\[line_1(3);line_2(3)];
        Intersect = [Intersect; new_intersect1'; new_intersect2'];
        
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