%% Intersection
% Intersection uses vector operations to determine the intersection point
% of two line segments. It was written with the assumption that the
% likelihood of the two segments NOT intersecting is very high, so it runs
% some computationally simple exclusionary tests first

% A and B are each defined by two vectors, their starting point and their direction.
% [Aposx Aposy; Adirx Adiry] ... same for B
% out is false if no intersection, [Ix Iy] if yes
% segments are given by their starting position vector and length/direction
% vector

% This is for 2D vectors only in the xy plane

function out = intersection(A,B)

%% The Far Away Exclusion
% first perform a simple test to exclude segments where both x or y are
% greater or less than both x or y of the other:

Axrange = [A(1,1) A(1,1)+A(2,1)];
Ayrange = [A(1,2) A(1,2)+A(2,2)];
Bxrange = [B(1,1) B(1,1)+B(2,1)];
Byrange = [B(1,2) B(1,2)+B(2,2)];

if outside(Bxrange,Axrange) || outside(Byrange,Ayrange)
    out = 0;
    disp('Far Away')
    return
end

%% Determining the Point

Apos = [A(1,:), 0]; %need to work in 3d for the cross func.
Adir = [A(2,:), 0];
Bpos = [B(1,:), 0];
Bdir = [B(2,:), 0];
posdiff = Bpos-Apos; %difference in positions
AxB = cross(Adir,Bdir); %the cross product of the direction vectors
poxdir = cross(posdiff,Adir); %posdiff x A direction to test colinearity

%% Check for Overlapping colinear lines

if AxB == 0 & poxdir == 0
    if dot(posdiff,Adir)>0 && dot(posdiff,Adir)<dot(Adir,Adir)
        disp('overlapping, A first')
        out = Bpos(1,1:2);
        return
    elseif dot(-posdiff,Bdir)>0 && dot(-posdiff,Bdir)<dot(Bdir,Bdir)
        disp('overlapping, B first')
        out = Apos(1,1:2);
        return
    else
        disp('colinear but disjoint')
        out = 0;
        return
    end
end

%% Determine point of intersection

t = cross(posdiff,Bdir)/AxB;
if t>1 || t<0
    out = 0;
    return
end

u = cross(posdiff,Adir)/AxB;
if u>1 || u<0
    out = 0;
    return
end

out = Bpos(1,1:2) + u*Bdir(1,1:2);

end

%% B Outside A?

% given two vectors that are the endpoints of ranges of values, is one
% fully outside the other? if yes, true, if no, false
% out is 0 or 1, default is 1

function out = outside(B,A)

out = 0;

if B(1) > max(A) && B(2) > max(A) || B(1) < min(A) && B(2) < min(A)
    out = 1;
    return
end
end
