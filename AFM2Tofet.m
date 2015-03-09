function [XYZ Edge] = AFM2Tofet(Imgpath)

disp('Segmenting')
[SEGfile,ORIENT,MAXCONF] = Full_Seg2(Imgpath);

disp('Binning')
bins=4;
[Binfile, MS] = Bin_Angles(SEGfile,bins);

disp('Molecularizing')
load(Binfile);
size = 400;
MS = MS(1:size,1:size,:);
pw = 2000/2160; ph = pw;
lsx = 0.37; lsy = 0.37; lsa = 0.38;

[XYZ, Edge] = Populate_Lattice(MS,pw,ph,lsx,lsy,lsa);

Latfile = [Imgpath(1:end-4), '_Lat'];

save(Latfile,'XYZ','Edge','-v7.3')

ImSizenm = size*(2000/2160);

[TofXY TofEdge] = AddEnergies(XYZ,Edge,Imgpath,ImSizenm);


end