function [TofXY TofEdge] = Bin2Tofet(Binfile)

load(Binfile);
size = 400;
MS = MS(1:size,1:size,:);
pw = 2000/2160; ph = pw;
lsx = 0.37; lsy = 0.37; lsa = 0.38;

[XYZ, Edge] = Populate_Lattice(MS,pw,ph,lsx,lsy,lsa);

Latfile = [Binfile(1:end-4), '_Lat'];

save(Latfile,'XYZ','Edge','-v7.3')

ImSizenm = size*(2000/2160);

[TofXY TofEdge] = AddEnergies(XYZ,Edge,ImSizenm);

end