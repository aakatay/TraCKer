% recursively finds the DS data folders and runs the DSrotate 
goRoot; % goes to the root folder
cd \MATLAB\TraCKer\
cd('data-biChang')
clear;
isClearfPath = 1;
isReConfig = 0;
dFolder = dir('data/DSdata/'); % 2013*  (dataFolder)

%fpath = rdir('/Volumes/My Passport/MATLAB/results from Bi-Chang scope/**/*DS','isdir==1');
if ismac
    fpath = rdir('/Volumes/My Passport/MATLAB/results from Wes scope/**/*DS','isdir==1');
    fRotated3D = '/Rotated3D/';
else %pc
    fpath = rdir('data\DSdata\**\*DS','isdir==1');
    fRotated3D = '\Rotated3D\'; 
end
if exist('fpathInd.mat')
    load('fpathInd.mat');
    fpath=fpath([fpathInd]);
    dFolder = dFolder([fpathInd]+2);
    jj=0;
else
    jj=4;
end
dF=0;
if ispc, dF = 4; end; % dummy folders
dF = 4;

if size(fpath,1) ~= size(dFolder,1)-dF
    display('WARNING: one of the folders has more than one data sets (DS)');
end    
CD = cd; % current directory

if ~exist(strcat(CD,'/Rotated3D'))
    mkdir(CD,'Rotated3D');
end

for i = 1 : size(fpath,1)
    cd(fpath(i).name)
    display(sprintf('---> %s \n', fpath(i).name));
    if 7~=exist(strcat(CD,fRotated3D,dFolder(i+jj).name))
        mkdir(strcat(CD,fRotated3D,dFolder(i+jj).name));
    end
    DSrotate;
    cd(CD);

end

%save('fpath.mat',fpath)
