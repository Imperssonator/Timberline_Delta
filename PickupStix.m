%% Pickup Stix
% This function will generate a structure array that will index all of the
% conducting structures in a given square, their origin, direction,
% and connections to other structures
% there will be two simultaneous data structures built: one that stores
% data on the fibers themselves, and the other that converts this data to
% nodes and edges

function out = PickupStix()

%% Nucleation
% generate the nucleii and store them in the fiber structures

length = 100; width = 100;                  %these values should correspond roughly to nanometers

F = struct;                                 % F is for Fibers. Each fiber is an element of a structure array
F = add_nucleii(F,length,width);
%disp(F)

%% Orientation
% pick an angle for growth for each nucleii

F = growth_directions(F);
%disp(F)

%% length
% pick a length from a distribution

F = grow_stacks(F);
disp(F)

end


%% Add Nucleii
% given a list of nucleus points and a fiber array F, add the nucleii as
% new fibers

function out = add_nucleii(F,l,w)

startpos = length(F);
N = nucleate(l,w);

for i = 1:length(N)
    F(i+startpos).nukept = N(i,:);
    F(i+startpos).type = 'stack';
end
out = F;
end



%% Nucleation function
% this will be random but I have allowed for any type of function to
% generate nucleii locations

function out = nucleate(length,width)

Area_density = 100;             % how many fibrils per unit area
out = rand(Area_density,2);
out = out * [length 0; 0 width];

end

%% Orient Growth
% for all 'stacks' in F, pick an angle from whatever distribution defines
% the fibers' orientation

function out = growth_directions(F)

for i = 1:length(F)
    if strcmp(F(i).type,'stack')
        F(i).angle = pick_angle();
    end
end
out = F;
end

%% Pick Angle
% no inputs, just a function that generates angles between 0 and 2pi from
% some distribution that I can specify within this function

function out = pick_angle()

ang = rand;
out = 2*pi*ang;

end

%% Grow Stacks
% given a structure, give a length to all the stacks that have nucleii and
% angles

function out = grow_stacks(F)

for i = 1:length(F)
    if strcmp(F(i).type,'stack') && getfield(F(i).nukept)~=[] && getfield(F(i).angle)~=[]
        F(i).length = pick_length();
    end
end
end

function out = pick_length()

len = randn;
avglen = 50;        % the average stack length in nm
stddev = 5;         % standard deviation of lengths
out = len*stddev+avglen;

end
       












