function MS = Bin_Angles(Segpath,bins)
%% Bin Angles
% will create an (image size) x # of microstates 3D binary matrix
% Each matrix will have 1's where that microstate exists

load(Segpath);

[m,n] = size(ORIENT);

MS = zeros(m,n,bins+1);

for i = 1:m
    for j = 1:n
        bin = ceil(ORIENT(i,j)/(180/bins))+bins/2;
        angle = ceil(ORIENT(i,j)/(180/bins))*(180/bins)-180/bins/2;
        if bin<1
            bin = bins+1;
        end
        MS(i,j,bin) = angle;
    end
end

MS = MS.*repmat(MAXCONF~=0,[1 1 bins+1]);
MS(:,:,bins+1) = MAXCONF==0;

end
    