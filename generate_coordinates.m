function [coord, bonds, type] = generate_coordinates(micstruct,pw,ph)
clc
% This function generates the xy coordinates of all the monomers present,
% list of which monomers are bonded to which and how
% The grid is taken to be in the 4th quadrant (origin top left)
% Assumes that the distance between two grid points is equal to the
% distance between two monomers (taken as 4 here)

[m, n] = size(micstruct);
coord = [];
bonds = [];
for i = 1:m
    for j = 1:n
        gridpt = [4*j 4*-i];                % gives coordinates in space based on position in matrix
        if micstruct(i,j) ~= 0
            coord = [coord; gridpt];        
        end
    end
end
% Calculating the distance between the current monomer and every other
% monomer before this
for i = 2:length(coord)
    d = [];
    for j = 1:i-1
        d(j) = sqrt(((coord(i,1) - coord(j,1))^2) + ((coord(i,2) - coord(j,2))^2)); 
    end
    if min(d) == 4                                  % Ensures the closest monomer present is bonded to it
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
for i = 1:length(bonds)
    if coord((bonds(i,1)),2) == coord((bonds(i,2)),2)
        type(i) = 1;
    else
        type(i) = 2;
    end
end
end 