function [pixelCoord, xyFiberSpace, monList, bonds, type] = generate_coordinates(micstruct,pw,ph)
clc
% This function generates the xy coordinates of all the monomers present,
% list of which monomers are bonded to which and how
% The grid is taken to be in the 4th quadrant (origin top left)

%% Hard Coded Parameters
% All parameters that are variable and hard-coded found here

piStackDistance = 0.37;     % Distance between two pi stacks in units of real space (3.7 A)
monLength = 0.37;           % Length along a chain that a single monomer occupies (also 3.7 A)

%% Variable Initializations
[m, n] = size(micstruct);
coord = [];
pixelCoord = [];            % Pixel Coord will store both the x and y cooridnates in real space and the indices of the pixel in the original image
xyFiberSpace = [];                            
% Fiber space is where the axis has been rotated by the angle in which the pi-stack is growing
% xyfiberspace gives the x & y coordinated in the rotated plane
bonds = [];
monList = [];               % Monomer List will store the real x-y coordinates of each monomer
type = [];

%% Conversion to Fiber Space relative to image origin
for i = 1:m
    for j = 1:n
        xyReal = [pw*j ph*-i];                % gives coordinates in space based on position in matrix
        if micstruct(i,j) ~= -180
            theta = micstruct(i,j);           % gives the angle with which axis should be rotated
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
stackLength = floor(L/piStackDistance)+1;                   % Number of chains in the stack rounded down

%% Monomer Placement

StartPixReal = pixelCoord(StartPix,1:2);
monList = [monList; [StartPixReal]];
ChainCenter = StartPixReal;
for i = 1:stackLength
    inFiber = 0;
    dir = 1;
    xCheck = 0;
    while not(inFiber)
        if IsInFiber(ChainCenter+[0 xCheck],rad,micstruct,pw,ph)
            ChainCenter = ChainCenter+[0 xCheck];
            monList = [monList; ChainCenter];
            break
        end
        if sign(xCheck)==1
            xCheck = xCheck*-1;
        else
            xCheck = xCheck*-1 + monLength;
        end
    end
    
    for dir = [1 -1];
        
        firstPlace = 1;
        inFiber = 1;
        
        while inFiber
            
            if firstPlace == 1
                prevPt = ChainCenter;
                firstPlace = 0;
            else
                prevPt = monList(end,:);
            end
            currPt = prevPt + [0 dir*monLength];
            if not(IsInFiber(currPt,rad,micstruct,pw,ph))
                inFiber = 0;
            else
                monList = [monList; currPt];
            end
        end
        
    end
    
    ChainCenter = ChainCenter + [piStackDistance 0];
    
end
    
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

function out = IsInFiber(FSpt,rad,micstruct,pw,ph)

RSpt = fs2real(FSpt,rad);
pixel_space = RSpt./[pw ph];
pixel_ind = ceil(pixel_space);
if micstruct(pixel_ind) == -180
    out = 0;
else
    out = 1;
end
end

function xyreal = fs2real(FSpt,rad)

xyreal = [0 0];
xyreal(1) = FSpt(1)*cos(rad) + FSpt(2)*sin(rad);
xyreal(2) = FSpt(1)*sin(rad) - FSpt(2)*cos(rad);

end



% function xyfs = real2fs(gridpt, rad)
% 
% out1 = [(gridpt(1)*cos(rad)+gridpt(2)*sin(rad)) (gridpt(2)*cos(rad)-gridpt(1)*sin(rad))]
% out2 = [dot(gridpt,[cos(rad) sin(rad)]), dot(gridpt,[cos(rad-pi/2) sin(rad-pi/2)])]
% 
% xyfs = 0;
% 
% end
% 


