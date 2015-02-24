load('MS_Example.mat');
MS = MS(1:200,1:200,:);
pw = 2000/2160; ph = pw;
lsx = 0.37; lsy = lsx;

[MonList, Edges] = Lattice_Fill(MS,pw,ph,lsx,lsy);

figure
spy(MS(:,:,1))
figure
plot(MonList(:,1),MonList(:,2),'.b')
