clear all
clc

Binfile = 'Tofet Files/15%_2min_Bin.mat';
load(Binfile);
size = 400;
MS = MS(1:size,1:size,:);
pw = 2000/2160; ph = pw;
lsx = 0.37; lsy = 0.37; lsa = 0.38;

[XYZ, Edge] = Populate_Lattice(MS,pw,ph,lsx,lsy,lsa);

% figure
% spy(MS(:,:,1))
% figure
% plot(MonList(:,1),MonList(:,2),'.b')
