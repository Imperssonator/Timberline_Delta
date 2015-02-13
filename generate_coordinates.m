function [XY, bonds, type] = generate_coordinates(MS,pw,ph)
clc
% This function generates the xy coordinates of all the monomers present,
% list of which monomers are bonded to which and how
% The grid is taken to be in the 4th quadrant (origin top left)

%% Hard Coded Parameters
% All parameters that are variable and hard-coded found here

piStackDistance = 0.37;     % Distance between two pi stacks in units of real space (3.7 A)
monLength = 0.37;           % Length along a chain that a single monomer occupies (also 3.7 A)

%% Variable Initializations
[m, n] = size(MS);
coord = [];
pixelCoord = [];            % Pixel Coord will store both the x and y cooridnates in real space and the indices of the pixel in the original image
xyFiberSpace = [];                            
% Fiber space is where the axis has been rotated by the angle in which the pi-stack is growing
% xyfiberspace gives the x & y coordinated in the rotated plane
bonds = [];
MonList = [];               % Monomer List will store the real x-y coordinates of each monomer
type = [];
ChainList = [];             % keep track of the length of each chain we form

%% Conversion to Fiber Space relative to image origin
for i = 1:m
    for j = 1:n
        xyReal = [pw*j ph*-i];                % gives coordinates in space based on position in matrix
        if MS(i,j) ~= -180
            theta = MS(i,j);           % gives the angle with which axis should be rotated
                                              % needs to be modified when we're dealing with more
                                              % than one angle
            rad = (theta/180)*pi;             % assuming the given angle is in degrees
            pixelCoord = [pixelCoord; [xyReal i j]];
            xyFiberSpace = [xyFiberSpace; [(xyReal(1)*cos(rad)+xyReal(2)*sin(rad)) (xyReal(2)*cos(rad)-xyReal(1)*sin(rad))]]; % basically the same as taking a dot product
            % xfs = xcos(theta) + ysin(theta)
            % yfs = ycos(theta) - xsin(theta)
        end
    end
end

[StackStart, StartPix] = min(xyFiberSpace(:,1));            % Return the Fiber space image origin coordinate of the starting pixel and which pixel it is
[StackEnd, EndPix] = max(xyFiberSpace(:,1));                % same for end
L = StackEnd-StackStart;                                    % Length of the pi stack
stackLength = floor(L/piStackDistance)+1;                   % Number of chains in the stack rounded down, +1 for the beginning

%% Monomer Placement

ChainCenter = xyFiberSpace(StartPix,:);                     % Initialize the Chain center as the fiber space coordinate at the start of the stack

for i = 1:stackLength
    
    MonList = [MonList; ChainCenter];
    inFiber = 1;
    prevPt = ChainCenter;                                   % previous point was the chain center
    MonListStart = size(MonList,1);                         % this keeps track of where in the MonList we started the current chain we're growing
    
    %% Add up
    while inFiber
        currPt = prevPt + [0 monLength];
        if not(IsInFiber(currPt,rad,MS,pw,ph))
            inFiber = 0;
        else
            MonList = [MonList; currPt];
            prevPt = currPt;
        end
    end
    
    %% Add down
    prevPt = ChainCenter;
    inFiber = 1;
    while inFiber
        currPt = prevPt + [0 -monLength];
        if not(IsInFiber(currPt,rad,MS,pw,ph))
            inFiber = 0;
        else
            MonList = [MonList; currPt];
            prevPt = currPt;
        end
    end
    
    %% Find new center for next chain
    ChainLen = size(MonList,1)-MonListStart;                                % Figure out how long the chain we just made is
    ChainList = [ChainList; ChainLen];                                      % Add it to the chain list, which keeps track of the chain lengths (for use in MW Distribution later...)
    for i = MonListStart:size(MonList,1)                                    
        if IsInFiber(MonList(i,:)+[piStackDistance 0],rad,MS,pw,ph)         % If the position one pi stack distance away from any monomer on the last chain we added is in the crystal...
            ChainCenter = MonList(i,:)+[piStackDistance 0];                 % just start the next chain there
            break
        elseif i == size(MonList,1)
            inFiber = 0;                                                    % If we check the whole chain and there's nowhere to go, then we're at the end of the fiber
        end
    end
end

%% Convert MonList back to real space
XY = zeros(size(MonList));
for i = 1:size(MonList,1)
    XY(i,:) = fs2real(MonList(i,:),rad);
end

%% Extra code for forming bonds later

%     
%     
%         % Go down
%         prevPt
%     
% 
% chainLength = [];
% store = xyFiberSpace;
% for i = 1:length(xyFiberSpace)
% 
%     for j = 1:length(store)
%         if xyFiberSpace(i,1) == store(j,1)
%             monList = [monList; store(j,:)];         % makes a list of the monomers located at the same length point
%             store(j,:) = [];
%         end
%     end
%     chainStart = monList(max(monList(:,2)),:);
%     chainEnd = monList(min(monList(:,2)),:);
%     chainLength = [chainLength; [((chainStart(2)-chainEnd(2))/monLength)+1  xyFiberSpace(i,1)]];
%     % [number of monomers in the chain, length point]
%     
% end    
% 
% % Calculating the distance between the current monomer and every other
% % monomer before this
% 
% for i = 2:length(pixelCoord)
%     d = [];
%     for j = 1:i-1
%         d(j) = sqrt(((pixelCoord(i,1) - pixelCoord(j,1))^2) + ((pixelCoord(i,2) - pixelCoord(j,2))^2)); 
%     end
%     if min(d) == monLength                          % Ensures the closest monomer present is bonded to it
%         for k = 1:length(d)
%             if d(k) == min(d)                       % Detects bonded monomers
%                 bonds = [bonds;[i k]];              % Adds to list of monomer pairs that are bonded
%             end
%         end
%     end
% end
% % Array that tells how each monomer pair in bonds is bonded.
% % 1- backbone bonding, 2- pi stack
% % THIS WORKS ONLY WHEN THE ANGLE IS 90
% 
% for i = 1:length(bonds)
%     if pixelCoord((bonds(i,1)),2) == pixelCoord((bonds(i,2)),2)
%         type(i) = 2;
%     else
%         type(i) = 1;
%     end
% end

end 

function out = IsInFiber(FSpt,rad,MS,pw,ph)
%% IsInFiber
% Checks if a point in fiber space is actually within the pixels of a fiber
% in real space

RSpt = fs2real(FSpt,rad);
% disp(RSpt)
pixel_space = fliplr(RSpt./[pw -ph]);   % 1st dimension of pixel space is the negative y-axis of real space...
pixel_ind = ceil(pixel_space);
[m,n] = size(MS);

if pixel_ind(1)<1 || pixel_ind(2)<1 || pixel_ind(1)>m || pixel_ind(2)>n
    out = 0;
elseif MS(pixel_ind(1), pixel_ind(2)) == -180
    out = 0;
else
    out = 1;
end
end

function xyreal = fs2real(FSpt,rad)
%% fs2real
% converts fiber space coordinates to real coordinates... rotation

xyreal = [0 0];
xyreal(1) = FSpt(1)*cos(rad) - FSpt(2)*sin(rad);
xyreal(2) = FSpt(1)*sin(rad) + FSpt(2)*cos(rad);

end

