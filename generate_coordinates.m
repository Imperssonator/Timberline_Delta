function [coord, bonds, type] = generate_coordinates(micstruct,pw,ph)
clc
% This function generates the xy coordinates of all the monomers present,
% list of which monomers are bonded to which and how
% The grid is taken to be in the 4th quadrant (origin top left)

[m, n] = size(micstruct);
coord = [];
xyfiberspace = [];                            
% Fiber space is where the axis has been rotated by the angle in which the pi-stack is growing
% xyfiberspace gives the x & y coordinated in the rotated plane
bonds = [];
for i = 1:m
    for j = 1:n
        gridpt = [pw*j ph*-i];                % gives coordinates in space based on position in matrix
        if micstruct(i,j) ~= -180
            theta = micstruct(i,j);           % gives the angle with which axis should be rotated
                                              % needs to be modified when we're dealing with more
                                              % than one angle
            rad = (theta/180)*pi;             % assuming the given angle is in degrees
            coord = [coord; gridpt];
            xyfiberspace = [xyfiberspace; [(gridpt(1)*cos(rad)+gridpt(2)*sin(rad)) (gridpt(2)*cos(rad)-gridpt(1)*sin(rad))]];
            % xfs = xcos(theta) + ysin(theta)
            % yfs = ycos(theta) - xsin(theta)
        end
    end
end
L = max(xyfiberspace(:,1)) - min(xyfiberspace(:,1));        % Length of the pi stack
piStackDistance = 2;
stackLength = round(L/piStackDistance);                     % Number of chains in the stack
monList = [];
chainLength = [];
monLength = 2;
for i = 1:length(xyfiberspace)
    for j = 1:length(xyfiberspace)
        if xyfiberspace(i,1) == xyfiberspace(j,1)
            monList = [monList; xyfiberspace(j,:)];         % makes a list of the monomers located at the same length point
        end
    end
    chainLength = [chainLength; [(max(monList(:,2))-min(monList(:,2)))/monLength  i]];
    % Calculates length of chanin at each length point along the pi stack
end    

% Calculating the distance between the current monomer and every other
% monomer before this

for i = 2:length(coord)
    d = [];
    for j = 1:i-1
        d(j) = sqrt(((coord(i,1) - coord(j,1))^2) + ((coord(i,2) - coord(j,2))^2)); 
    end
    if min(d) == monLength                          % Ensures the closest monomer present is bonded to it
        for k = 1:length(d)
            if d(k) == min(d)                       % Detects bonded monomers
                bonds = [bonds;[i k]];              % Adds to list of monomer pairs that are bonded
            end
        end
    end
end
% Array that tells how each monomer pair in bonds is bonded.
% 1- backbone bonding, 2- pi stack
% THIS WORKS ONLY WHEN THE ANGLE IS 90
type = [];
for i = 1:length(bonds)
    if coord((bonds(i,1)),2) == coord((bonds(i,2)),2)
        type(i) = 2;
    else
        type(i) = 1;
    end
end
end 