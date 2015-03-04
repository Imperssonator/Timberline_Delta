NNBoxWidth = 0.8;  % Width of box to search for nearest neighbors.. >2*lattice spacing
MaxEdgeDist = 0.4; % maximum edge length that we will allow
MinEdgeDist = 0.3;

% disp(size(XYZ))
[XYZ, Edge] = FindNN(XYZ,Edge,NNBoxWidth,MaxEdgeDist,MinEdgeDist);