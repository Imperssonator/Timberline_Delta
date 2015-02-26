function [OuterMonList, Edges] = Lattice_Fill(MS,pw,ph,lsx,lsy)

%% Lattice Fill
% MS is a binned microstructure where each pixel in a bin represents its
% angle off the x-axis

% pw = pixel width
% ph = pixel height
% lsx = lattice spacing x
% lsy = lattice spacing y

MonList = [];
[h, w, bins] = size(MS);
h = h*ph;
w = w*pw;

%% Generate Monomer Grid


    
OuterMonList = []; % OuterMonList builds up the MonLists from each bin
for i = 1:1 % Right now just working on bin 1
    disp(i)
    Img_Corners = [0 0;...
        w 0;...
        0 h;...
        w h]; % Edges of the frame marked
    MSi = MS(:,:,i);
    Bin_Angle = MSi(find(MSi,1))/180*pi; % Finds the first non-zero value in the bin
    if isempty(Bin_Angle)
        break
    end
    Img_Corners_Rot = Rotate_List(Img_Corners,Bin_Angle); % Rotate the frame by bin angle
    
    [max_y ,Top_Corner] = max(Img_Corners_Rot(:,2));
    [min_y ,Bottom_Corner] = min(Img_Corners_Rot(:,2));
    [max_x ,Right_Corner] = max(Img_Corners_Rot(:,1));
    [min_x ,Left_Corner] = min(Img_Corners_Rot(:,1));
    
    % Equations for the lines of the box
    X_bottomright = [Img_Corners_Rot(Bottom_Corner,1) Img_Corners_Rot(Right_Corner,1)];
    Y_bottomright = [Img_Corners_Rot(Bottom_Corner,2) Img_Corners_Rot(Right_Corner,2)];
    line_bottomright = polyfit(X_bottomright, Y_bottomright,1); % Fit the two vertices of the frame to give the edge
    % y = ax + b, line_bottomright = [a b]
    
    X_topright = [Img_Corners_Rot(Top_Corner,1) Img_Corners_Rot(Right_Corner,1)];
    Y_topright = [Img_Corners_Rot(Top_Corner,2) Img_Corners_Rot(Right_Corner,2)];
    line_topright = polyfit(X_topright, Y_topright,1);
    
    X_topleft = [Img_Corners_Rot(Top_Corner,1) Img_Corners_Rot(Left_Corner,1)];
    Y_topleft = [Img_Corners_Rot(Top_Corner,2) Img_Corners_Rot(Left_Corner,2)];
    line_topleft = polyfit(X_topleft, Y_topleft,1);
    
    X_bottomleft = [Img_Corners_Rot(Bottom_Corner,1) Img_Corners_Rot(Left_Corner,1)];
    Y_bottomleft = [Img_Corners_Rot(Bottom_Corner,2) Img_Corners_Rot(Left_Corner,2)];
    line_bottomleft = polyfit(X_bottomleft, Y_bottomleft,1);
   
    Intersect = [Img_Corners_Rot(Bottom_Corner,:)];
    % gives the list of points on the edges of the frame
    crystalxy = zeros(ceil((h/lsy)*(w/lsx)),2); % gives the list of points that fills the lattice in crystal space
%     whos crystalxy
    crystalxy(1,:) = Intersect./[lsx lsy];
    county = 1;
    disp(min_y)
    disp(max_y)
    
    for y = min_y+lsy:lsy:max_y
        %         disp(y)
        line_1 = [0 1 y]; % ax + by = c for the y line
        
        % Right intersection
        if y<=Img_Corners_Rot(Right_Corner,2)
            line_2 = [line_bottomright(1) -1 -1*line_bottomright(2)]; % line_2 = [a b c] for the edge of the frame
        else
            line_2 = [line_topright(1) -1 -1*line_topright(2)];
        end
        
        new_intersect1 = [line_1(1:2); line_2(1:2)]\[line_1(3);line_2(3)]; % POI of line1 and line 2
        crystal_edge1 = new_intersect1'./[lsx lsy]; 
        
        % Left intersection
        if y <= Img_Corners_Rot(Left_Corner,2)
            line_2 = [line_bottomleft(1) -1 -1*line_bottomleft(2)];
        else
            line_2 = [line_topleft(1) -1 -1*line_topleft(2)];
        end
        
        new_intersect2 = [line_1(1:2); line_2(1:2)]\[line_1(3);line_2(3)];
        crystal_edge2 = new_intersect2'./[lsx lsy];
        
        Intersect = [Intersect; new_intersect1'; new_intersect2'];
        
        for x = ceil(crystal_edge2(1)):floor(crystal_edge1(1))
            county = county+1;
            crystalxy(county,:) = [x y]; % Filling the box
        end
        
        
    end
    
    crystalxy = crystalxy(1:county,:); % remove tail-end zeros
    xy_Rot = ones(size(crystalxy)); % Lattice points in rotated space (nm)
    disp('Filled in box')
    for j = 1:length(xy_Rot)
        xy_Rot(j,:) = crystalxy(j,:).*[lsx lsy]; % Scaling up to nm
    end
    Lattice_list = Straighten_List(xy_Rot, Bin_Angle); % Brings back frame to real space
    disp('Rotated box back')
    MonList = zeros(size(Lattice_list));
    countm = 0;
    
    for j = 1:length(Lattice_list)
        %         tic
        if IsInFiber(Lattice_list(j,:),MSi,pw,ph)
            countm = countm+1;
            MonList(countm,:) = Lattice_list(j,:); % Checks if each lattice point has a monomer
        end
        %         checker = toc
    end
    MonList = MonList(1:countm,:); % Remove tail-end zeros from preallocation
    Edges = Straighten_List(Intersect, Bin_Angle); % unnecessary but need an output for 'Edges' 
    OuterMonList = [OuterMonList; MonList];
    whos OuterMonList  % whos checks the size and type of any variable
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

function Straightened = Straighten_List(Rot_List,angle)

Straight = [cos(angle) sin(angle); % This should rotate the frame back by the same amount
            -sin(angle) cos(angle)];
Straightened = zeros(size(Rot_List));

for i = 1:size(Rot_List,1)
    Straightened(i,:) = Straight*Rot_List(i,:)';
end

end

function out = IsInFiber(RSpt,MS,pw,ph)
%% IsInFiber
% Checks if a point in fiber space is actually within the pixels of a fiber
% in real space

pixel_space = fliplr(RSpt./[pw ph]);   % 1st dimension of pixel space is the negative y-axis of real space...
pixel_ind = ceil(pixel_space);
[m,n] = size(MS);

if pixel_ind(1)<1 || pixel_ind(2)<1 || pixel_ind(1)>m || pixel_ind(2)>n
    out = 0;
elseif MS(pixel_ind(1), pixel_ind(2)) == 0 % If this real space point is not in this bin, remove from list
    out = 0;
else
    out = 1;
end
end