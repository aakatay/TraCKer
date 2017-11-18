function DSrotateCore(fNameRd,fNameWr,fPath,CD,dFolder,imgSize,rotAngle,zRange12,isReConfig)
% rotate the 3D images
    
    %Pre-size the movie matrix
    imgSize = [imgSize(1) imgSize(2) imgSize(3)];
    %imgSize = [2,21,10];

    for z=1:imgSize(3) % # of frames
        img3D(:,:,z) = imread(char(strcat(fPath,fNameRd)),z); 
        %img3D = flipdim(img3D,2);
    end

    img3D = permute(img3D,[2 3 1 ]);
    img3D = imrotate(img3D,rotAngle,'bicubic');
    img3D = permute(img3D,[3 1 2 ]);

    %% dorsal - ventral separation
    Z=size(img3D,3);
    %load('DV.mat');   %DV = DV + ( Z - zRange12(2)-1 ); 
    DV = 100;
    %% calculate lateral tilt angle and update DV (first image only)
    if isempty(strfind(char(fNameWr),'_001.')) && isempty(strfind(char(fNameWr),'_01.')) % first image
        ;
    else
        %imagesc(maxPventral);
        %figure
        if ~exist('lateralTilt.mat') || isReConfig
            Np = 5; % number of points for each frame
            in = 1;
            %construct image base stack
            fr1 = 3; % # of frames before DV
            fr2 = 15; % # of frames after DV
            imgBase = img3D(:,:,DV-fr1 : DV+fr2);
            q = 0; xyS =[0 0]; zS =0; updDV =0; isD = 0; isV =0; wSt = 0;
            while q == 0
                [mx, mxind] = max(imgBase(:));
                imgBase = imgBase/mx*128;
                [X_, Y_, zIx] = ind2sub(size(imgBase),mxind);
                nS = 1;
                % graphical input of the spot positions
                while (1)
                    imgSec = imgBase(:,:,zIx); 
                    image(imgSec); axis equal;
                    h = figure(1);
                    btn = 0;
                    while btn == 0
                        btn = waitforbuttonpress;
                        k = get(h,'CurrentCharacter');
                    end
                    switch lower(k)
                        case 'a'
                            zIx = zIx - 1;
                        case 'd'
                            zIx = zIx + 1;
                        case 's' % select spots
                            xyS(nS,:) = ginput;
                            zS(nS) = zIx;  % +DV-fr1-1
                            zIx = zIx + 1; % index for selected spots
                            nS = nS + 1;
                        case 'w'  % check again
                            xyz = floor([xyS'; zS]);                            
                            phiTheta = findThetaPhi(imgBase,xyz)
                            imgBase = correctLateralTilt(imgBase,phiTheta);  
                            zIx = DV;   
                            wSt = 1;                       
                            break
                        case 'q' % quit tilt correction
                            q = 1;
                            break
                        case 'e' % load the whole stack
                            if ~exist('phiTheta')
                                phiTheta = [0,0];
                                imgBase = img3D;
                            else
                                imgBase = correctLateralTilt(img3D,phiTheta); 
                            end
                            [mx, mxind] = max(imgBase(:));
                            [X_, Y_, zIx] = ind2sub(size(imgBase),mxind);
                            zIx = DV;
                            wSt = 1;
                            break
                        case 'v' % update DV value
                            updDV = 1;
                            DVupd = zIx;
                        case '1'
                            isD = 1;
                            crop1 = zIx; % dorsal
                            display(sprintf('dorsal crop selected: %d',crop1));
                        case '2'
                            isV = 1;
                            crop2 = zIx; % ventral
                            display(sprintf('ventral crop selected: %d',crop2));
                    end
                    if lower(k) == 's' && (xyS(nS) == 0 || xyS(nS) == 0)
                        zIx = zIx - 1;
                        nS = nS - 1;
                        q = 0;
                        display('wrong input detected, re-enter');
                    end
                    %% keep the z-scan in the valid range
                    if wSt == 1 %whole stack
                        if zIx > size(imgBase,3)
                            zIx = size(imgBase,3);
                        elseif zIx < 1
                            zIx = 1;
                        end
                    else
                        if zIx > fr1+fr2+1
                            zIx = fr1+fr2+1;
                        elseif zIx < 1
                            zIx = 1;
                        end
                    end
                    if isD && isV % define crop range
                        zRange12 = [crop1 crop2];
                        save('zRange12.mat','zRange12');
                        display('zRange updated');
                        isD=0; isV=0;
                    end
                end % while (1)
            end % while (quit)
            if updDV
                DV = DVupd
                display('DV updated')
                save('DV.mat','DV');
            end
            save('lateralTilt.mat','phiTheta')
        else
            load('lateralTilt.mat')
        end
    end
    load('lateralTilt.mat') 
    %% correct for lateral tilt
    if phiTheta(1) ~= 0 || phiTheta(2) ~= 0 
        img3D = correctLateralTilt(img3D,phiTheta);
    end
    

    % clip edges from right
    clipEdge = 6; % [px] clip images from the right 
    for i = 1:size(img3D,3)
        if clipEdge
            [I,ix1]=find(img3D(1,:,i)>0); 
            [I,ix2]=find(img3D(end,:,i)>0); 
            ix = max([ix1 ix2]);
            clipEdge2 = clipEdge+abs(max(ix2)-max(ix1));
            if isempty(clipEdge2)
                ;
            else
                img3D(:,ix+1-clipEdge2:max(ix),i)=0;
            end
        end
    end
    %% 3D
    if 7~=exist(strcat(CD,'/data/Rotated3D'))
        mkdir(CD,'/data/Rotated3D');
    end
    if 7~=exist(strcat(CD,'/data/Rotated3D/',dFolder))
        mkdir(strcat(CD,'/data/Rotated3D/'),dFolder);
    end
    if 7~=exist(strcat(CD,'/data/Rotated3D/',dFolder,'/3D'))
        mkdir(strcat(CD,'/data/Rotated3D/',dFolder),'3D');
    end
    fnameRtd = strcat(CD,'/data/Rotated3D/',dFolder,'/3D/',fNameWr);
    writeTiffStack(img3D,fnameRtd,zRange12,DV,3);
    %% dorsal-ventral 3D
    if 7~=exist(strcat(CD,'/data/Rotated3D/',dFolder,'/dorsal3D'))
        mkdir(strcat(CD,'/data/Rotated3D/',dFolder),'dorsal3D');
    end
    if 7~=exist(strcat(CD,'/data/Rotated3D/',dFolder,'/ventral3D'))
        mkdir(strcat(CD,'/data/Rotated3D/',dFolder),'ventral3D');
    end    
    fnameDorsal3D = strcat(CD,'/data/Rotated3D/',dFolder,'/dorsal3D/',fNameWr);
    %writeTiffStack(img3D,fnameDorsal3D,zRange12,DV,1);
    fnameVentral3D = strcat(CD,'/data/Rotated3D/',dFolder,'/ventral3D/',fNameWr);
    %writeTiffStack(img3D,fnameVentral3D,zRange12,DV,2);
    %% max-projection
    maxP        = max(img3D(:,:,zRange12(1):zRange12(2)),[],3);
    maxPdorsal  = max(img3D(:,:,zRange12(1):DV-1),[],3);
    [maxPventral MxZeInd]= max(img3D(:,:,DV:zRange12(2)),[],3);
    
    fnameMxProj         = strcat(CD,'/data/Rotated3D/',dFolder,'/3D/','maxProj-',fNameWr);
    fnameMxProjDorsal   = strcat(CD,'/data/Rotated3D/',dFolder,'/dorsal3D/','maxProj-',fNameWr);
    fnameMxProjVentral  = strcat(CD,'/data/Rotated3D/',dFolder,'/ventral3D/','maxProj-',fNameWr);
    imwrite(uint16(maxP),char(fnameMxProj),'tiff');
    imwrite(uint16(maxPdorsal),char(fnameMxProjDorsal),'tiff');
    imwrite(uint16(maxPventral),char(fnameMxProjVentral),'tiff');
    

    
end

%[maxPventral MxZeInd]= max(img3D(:,:,DV:DV+5),[],3);  figure(1); imagesc(maxPventral); figure(2); imagesc(MxZeInd);
%[maxPventral MxZeInd]= max(img3D(1:100,1:100,DV:DV+5),[],3); figure(1); imagesc(maxPventral); figure(2); imagesc(MxZeInd);
