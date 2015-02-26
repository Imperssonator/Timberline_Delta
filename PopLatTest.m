clear all
clc

load('MS_Example.mat');
size = 400;
MS = MS(1:size,1:size,:);
pw = 2000/2160; ph = pw;
lsx = 0.37; lsy = 0.7;

[XYZ, Edge] = Populate_Lattice(MS,pw,ph,lsx,lsy);

% figure
% spy(MS(:,:,1))
% figure
% plot(MonList(:,1),MonList(:,2),'.b')
