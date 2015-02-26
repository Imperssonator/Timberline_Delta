function [XYZ, Edge] = Populate_Lattice(MS,pw,ph,lsx,lsy)

%% Populate Lattice
% This function takes a binned 2-D microstructure MS, where each pixel in a
% bin is either '0' or the angle of its crystalline domain, and creates a
% list of molecules that completely fill all grains in that microstructure
% with the correct lattice orientations and spacings.
%
% pw = pixel width in nm
% ph = pixel height in nm
% lsx = lattice spacing perpendicular to the angle of the grain
% lsy = lattice spacing along the axis of the grain specified by its angle

% XYZ will be a list of molecule locations along with their site energies
% Edge will be a list of which molecules in XYZ are connected to each other
% for the purposes of charge transport, and the weight of that edge.

%% Important Hard Coded Parameters

Backbone = 1; % These will be used to denote the type of bond between each site
PiStack = 2;
Amor = 3;

%% Find Frame Edges

[Ih Iw] = size(MS(:,:,1)); % image height and width in pixels
nmh = Ih*ph; nmw = Iw*pw; % image height and width in nm
TopLeft = [0 0];           % The corners of our image in real space (nm)
TopRight = [nmw 0];
BottomLeft = [0 -nmh];
BottomRight = [nmw -nmh];

%% Initialize Variables
XYZ = [];
Edge = [];

%% Populate the molecules of each crystalline bin

bins = size(MS,3);

for bin = 1:bins-1
    MSbin = MS(:,:,bin); % Generate the microstructure for this bin
    BinAngle = MSbin(find(MSbin,1))/180*pi; % The bin angle is stored at every pixel, so just find the first one
    if isempty(BinAngle)
        break
    end
    
    %% Compute Scaling and Shifting Parameters
    if BinAngle<=0
        Diag = (nmw^2+nmh^2)^(1/2); % Image diagonal
        LatticeWidth = Diag*cos(BinAngle+pi/4); % How wide of a square at angle BinAngle do we need to fully circumscribe the image
        ShiftVector = nmw*sin(-BinAngle).*[cos(BinAngle+pi/2) sin(BinAngle+pi/2)];
    else
        Diag = (nmw^2+nmh^2)^(1/2);
        LatticeWidth = Diag*cos(pi/4-BinAngle); % How wide of a square at angle BinAngle do we need to fully circumscribe the image
        ShiftVector = -nmh*sin(BinAngle).*[cos(BinAngle) sin(BinAngle)];
    end
    
    %% Make a mesh grid
    disp('Making Grid')
    [FiberGrid ChainGrid] = ndgrid(0:lsy:LatticeWidth,0:-lsx:-LatticeWidth); % The "Y" axis in this case is the fiber axis and starts off as the real space x-axis
    
    disp('Populating Edges')
    [NumChains ChainLen] = size(FiberGrid);
    EdgeBin = zeros(2*NumChains*ChainLen, 3); % Initialize Edges so we can fill them in for each grain
    GridPts = [FiberGrid(:) ChainGrid(:)];
    NumPtsOrig = size(GridPts,1);
    %% Populate Edges
    EdgeCount = 1;
    for m = 0:ChainLen-2
        for c = 1:NumChains-1
            EdgeBin(EdgeCount,:) = [m*NumChains+c m*NumChains+c+1 PiStack]; % Fill in Pi-Stack bonds
            EdgeCount = EdgeCount+1;
            
            EdgeBin(EdgeCount,:) = [m*NumChains+c (m+1)*NumChains+c Backbone]; % Fill in Backbone bonds
            EdgeCount = EdgeCount+1;
        end
    end
    for m = 1:ChainLen-1
        EdgeBin(EdgeCount,:) = [m*NumChains (m+1)*NumChains Backbone]; % Last row of backbones
        EdgeCount = EdgeCount+1;
    end
    for c = 1:NumChains-1
        EdgeBin(EdgeCount,:) = [(ChainLen-1)*NumChains+c (ChainLen-1)*NumChains+c+1 PiStack];
        EdgeCount = EdgeCount+1;
    end
    
    EdgeBin(EdgeBin(:,1)<=0,:)=[]; % Clean it up
    
    %% Rotate and cut the lattice
    disp('Rotating Grid')
    RotGridPts = Rotate_List(GridPts,BinAngle);
    disp('Shifting Grid')
    ShiftedGP = RotGridPts + repmat(ShiftVector,[size(RotGridPts,1) 1]);
    ShiftedGP = [ShiftedGP (1:size(ShiftedGP,1))'];
    disp('Cutting Grid')
    ShiftedGP(ShiftedGP(:,1)>nmw,:)=[];
    ShiftedGP(ShiftedGP(:,1)<0,:)=[];
    ShiftedGP(ShiftedGP(:,2)>0,:)=[];
    ShiftedGP(ShiftedGP(:,2)<-nmh,:)=[];
    
    %% For every point on the lattice, keep it if it's in this bin's domain
    disp('Filling Domains')
    NumPtsCut = size(ShiftedGP,1);
    for i = 1:NumPtsCut
        CheckPt = ceil([ShiftedGP(i,2)/-ph ShiftedGP(i,1)/pw]);
        if MSbin(CheckPt(1),CheckPt(2)) == 0
            ShiftedGP(i,:) = [-1 0 ShiftedGP(i,3)];
        end
    end
    ShiftedGP(ShiftedGP(:,1)<0,:)=[];
    
    disp('Fixing Edges')
    NumEd = size(EdgeBin,1); % Original number of edges
    NumPtsNew = size(ShiftedGP,1); % Number of Points after cutting and filling
%     FixedEdge = zeros(size(EdgeBin)); % A list to rebuild the edges now that points are gone
%     EdCount = 1; % Keep track of how many have been added
    PtShift = zeros(NumPtsOrig,1);  % A handy list to show where original grid points ended up
    
    for i = 1:NumPtsNew             % i = # of new point
        % Take all the edges that involve that point, add them to the fixed
        % list, and remove them from the original list
        % Crucial point here: edges always connect from lower to higher
        op = ShiftedGP(i,3);        % original point
        PtShift(op,1) = i;          % store where the old point went
    end
    
    for j = 1:NumEd
        if PtShift(EdgeBin(j,1),1) == 0 || PtShift(EdgeBin(j,2),1) == 0 % If either original point didn't make the new list
            EdgeBin(j,:) = [-1 -1 -1];
        else
            EdgeBin(j,:) = [PtShift(EdgeBin(j,1),1) PtShift(EdgeBin(j,2),1) EdgeBin(j,3)];
        end
    end
    EdgeBin(EdgeBin(:,1)<0,:)=[];
%     
%     figure
%     plot(ShiftedGP(:,1),ShiftedGP(:,2),'.b')
    XYZ = [XYZ; ShiftedGP];
    Edge = [Edge; EdgeBin];
end

figure
plot(XYZ(:,1),XYZ(:,2),'.b')
axis equal

end

function Rotated = Rotate_List(PtList,angle)

ROT = [cos(angle) -sin(angle);...
       sin(angle) cos(angle)];
Rotated = zeros(size(PtList,1),2);
   for i = 1:size(PtList,1)
       Rotated(i,:) = ROT*PtList(i,:)';
   end
end

function out = IsInFiber(RSpt,MS,pw,ph)
%% IsInFiber
% Checks if a point in fiber space is actually within the pixels of a fiber
% in real space

pixel_space = fliplr(RSpt./[pw -ph]);   % 1st dimension of pixel space is the negative y-axis of real space...
pixel_ind = ceil(pixel_space);

if MS(pixel_ind(1), pixel_ind(2)) ~= 0 % If this real space point is not in this bin, remove from list
    out = 1;
else
    out = 0;
end
end

