function command = Bin2Tofet(Dir,ImName,StartPix,BoxSize)
% Dir is the directory where the AFM folders are located
% StartPix is a 1x2 vector of integers where the simulation box should
% start
% Size is the size of the box edges in pixels

% Dir should be 'Tofet Files/'
% ImName should be '15%_2min'

%% Generate File Locations
SubDir = [Dir ImName '/'];                  % Tofet Files/15%_2min/
Binfile = [SubDir ImName '_Bin'];           % Tofet Files/15%_2min/15%_2min_Bin
PixStr = mat2str(StartPix);                 % '[10 10]'
PixStr = PixStr(2:end-1);                   % '10 10'
SizeStr = mat2str(BoxSize);                    % '50'
SizeName = [ImName '_' PixStr '_' SizeStr]; % '15%_2min_10 10_50'
SizeDir = [SubDir SizeName '/']; % 'Tofet Files/15%_2min/15%_2min_10 10_50/'
mkdir(SizeDir)

%% Cut the original bin to the specified size

load(Binfile);

MS = MS(StartPix(1):StartPix(1)+BoxSize,StartPix(2):StartPix(2)+BoxSize,:);

%% Run Populate Lattice

pw = 2000/2160; ph = pw;
lsx = 0.37; lsy = 0.37; lsa = 0.38;

[XYZ, Edge] = Populate_Lattice(MS,pw,ph,lsx,lsy,lsa);

Latfile = [SizeDir SizeName '_Lat'];        % 'Tofet Files/15%_2min/15%_2min_10 10_50/15%_2min_10 10_50_Lat'

save(Latfile,'XYZ','Edge','-v7.3')

%% Convert the Lattice to .xyz and .edge files

ImSizenm = BoxSize*pw;

[XYZFile EdgeFile] = AddEnergiesTof(XYZ,Edge,SizeDir,SizeName,ImSizenm);

%% Run Tofet

tft = '''/Users/Imperssonator/CC/Tofet/examples/GSL_randomGenerator/tof/tft'' ';
AD = pwd;
sim = '''/Users/Imperssonator/CC/Tofet/examples/GSL_randomGenerator/AFMSims/P3HTtof.sim'' ';
xyz = ['''' AD '/' XYZFile ''' '];
edge = ['''' AD '/' EdgeFile ''' '];
out = ['''' AD '/' SizeDir SizeName '.out'''];

% 
% command = [tft sim xyz edge '> ' out];
% system(command)

end