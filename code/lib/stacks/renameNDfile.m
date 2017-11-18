% converts ND2 file to tiff file 
% names according to:
% filename001.ND2 --> 76x60xy_18ms_100P_300G_1x_003
% xy: crop region
% ms: exposure time
% P: power [percent]
% G: gain
% x: magnification

clear
fname = dir('*.nd2');
nF = numel(fname);

for i = 1:nF
    clear data;
    fn = fname(i).name;
    data = bfopen(fn);
    metadata = data{1, 2};
    A = data{1,1};
    metadataKeys = metadata.keySet().iterator();
    
    keyVal(1) = 100; % default power level
    for i=1:metadata.size()
        key = metadataKeys.nextElement();
        if strfind(key,'Global Z position for position, plane'), continue;end;
        if strfind(key,'Z position for position, plane'), continue;end;
        if strfind(key,'timestamp #'), continue;end;        
        key2 = key(8:end);

        keys = {'Global Line','Exposure','dMinPeriodDiff','dMaxPeriodDiff','dAvgPeriodDiff','GainMultiplier','PFS, offset','CurrentTemperature','right','left','top','bottom','top','Zoom'};
        keysNum = {'dMinPeriodDiff','dMaxPeriodDiff','dAvgPeriodDiff'}; % keys with numeric values
                    
        ix = find(strcmp(key2,keys) == 1);
        if strfind(key2,'Line')
            value = metadata.get(key);
            if strfind(value,'Active')
                keyVal(1)= str2num(value(1:end-8));
            end      
        elseif ~isempty(ix)
            if strcmp(key2,'Zoom') == 1
                temp = metadata.get(key);
                keyVal(ix) = str2num(temp(1:4));
            elseif find(strcmp(key2,keysNum) == 1)
                keyVal(ix) = metadata.get(key);
            else
                keyVal(ix) = str2num(metadata.get(key));
            end
        end
    end
    
%% generate filename [tif]
    sX = abs(keyVal(11)-keyVal(12)); % size x
    sY = abs(keyVal(9)-keyVal(10)); % size y
    fOut = sprintf('%ix%ixy_%ims_%iP_%iG_%ix_%s.tif',sX,sY,keyVal(2),keyVal(1),keyVal(6),keyVal(14),fn(end-6:end-4));

%% shape data
    nZ = 1; % # of z planes
    s1 = 1; % first frame
    s2 = nZ; % last frame
    h = waitbar(0,'reading data stack...');
    for i = 1 : size(A,1)/nZ
        for j = s1 : s2 % selected frames
            temp = cell2mat(A((i-1)*nZ+j,1));
            %Zstack(:,:,j) = temp(200:400,150:350); % CROP 
            Zstack(:,:,j-s1+1) = temp;
        end
        if nZ >1
            fn = sprintf('%s%02i.tif', filenameStack_out,i);
        end
        maxProj(:,:,i) = max(Zstack,[],3);
        waitbar(i / size(A,1)*nZ)
    end
    
%% write tiff file
    stackWrite(maxProj,fOut);
    close(h)
end