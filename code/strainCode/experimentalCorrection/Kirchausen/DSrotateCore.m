function DSrotateCore(fNameRd,fNameWr,imgSize,rotAngle,pxZvsX,isReConfig)
% rotate the 3D images
    
    for z=1:imgSize(3) % # of frames
        img3D_(:,:,z) = imread(fNameRd,z); 
        %img3D = flipdim(img3D,2);
    end
        
    img3D = permute(img3D_,[2 3 1 ]);
    %stackWrite(img3D,'img3D1.tif');
    if 1
        for i=1:size(img3D,3)
            x = 1:size(img3D,2); y = 1:size(img3D,1);
            xq = [1:(size(img3D,2)*pxZvsX)]/pxZvsX; yq =y;
            [Xq,Yq]=meshgrid(xq,yq);
            [X,Y]=meshgrid(x,y);
            img3Dintp(:,:,i) = interp2(X,Y,img3D(:,:,i),Xq,Yq);
        end
    end
    img3D = imrotate(img3Dintp,rotAngle,'bicubic');
    %stackWrite(img3D,'img3D2.tif');
    img3D = permute(img3D,[3 2 1 ]);
    %stackWrite(img3D,'img3D3.tif');
    
    %% dorsal - ventral separation
    Z=size(img3D,3);
    if exist('DV.mat'), load('DV.mat'); end;  
    %% calculate lateral tilt angle and update DV (first image only)
    if ~isempty(strfind(char(fNameWr),'_001.')) || ~isempty(strfind(char(fNameWr),'_01.')) % only for the first image
        if ~exist('lateralTilt.mat') || ~exist('DV.mat') || isReConfig % generate mat files:  DV lateralTilt
            Np = 5; % number of points for each frame
            in = 1;
            %construct image base stack
            if exist('DV') % DV loaded
                fr1 = DV(2)-3; % # of frames before Dv
                fr2 = DV(3)+10; % # of frames after dV
            else
                fr1 = 1;
                fr2 = Z;
            end
            description = sprintf('1-press e to load full stack. \n 2-correct for lateral tilt. \n 3-select D-V regions \n 4-quit');
            imgBase = img3D(:,:,fr1 : fr2);
            q = 0; xyS =[0 0]; zS =0; updDV =0; wSt = 0;
            h = figure(1);
            maximize
            gry = 0.8;
            yT = 140;
            hTextDirections = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,...
                'String',sprintf(' d: stack+1 \n a: stack-1 \n s: select spots \n w:reload \n e: load full stack \n REGISTER frame as: \n1: dorsal \n 2: DV (dorsal) \n 3: DV (ventral) \n 4: ventral '),...
                'Position',[0 100 yT 140],...
                'HorizontalAlignment','Left');  
            btn = 0;
            k=0;
            isSelSpots = 0;
            while q == 0 % quit is not pressed
                [mx, mxind] = max(imgBase(:));
                imgBase = imgBase/mx*128;
                [~, ~, zIx] = ind2sub(size(imgBase),mxind);
                nS = 1;
                if ~exist('DV'), DV(1)=1; DV(4)=Z; DV(2)=zIx; DV(3)=zIx+1; end;
                % graphical input of the spot positions
                if k == 'e', state = 'reloading'; else state = 'loading'; end;
                while (1)
                    if isSelSpots
                        xyS(nS,:) = ginput;
                        zS(nS) = zIx;  % +DV-fr1-1
                        zIx = zIx + 1;
                        nS = nS + 1; % index for selected spots
                        state = 'select spots to correct for lateral tilt';
                        if nS >1
                            description = sprintf('1-) browse to the stack layer where the spot has max intensity \n2-) click on the spot center and press ENTER \n3-) do the same some other spots in different regions');                
                        else % nS = 1
                            description = sprintf('select another spot or reload for tilt correction');
                        end
                        isSelSpots =0;
                    end
                    hTextStackNo = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,...
                        'String',sprintf(' stack# : %i ',zIx),...
                        'Position',[0 yT+150 140 30],...
                        'HorizontalAlignment','Left');  
                    hTextState = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,...
                        'String',sprintf(' state : %s ',state),...
                        'Position',[0 yT+200 140 30],...
                        'HorizontalAlignment','Left');  
                    hTextDescription = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,...
                        'String',sprintf(' description :\n%s ',description),...
                        'Position',[0 yT+250 140 150],...
                        'HorizontalAlignment','Left');
                    hTextDV = uicontrol('style','text','BackgroundColor',[1 1 1]*gry,...
                        'String',sprintf(' DORSAL-VENTRAL POSITIONS: \n  ie. (D >= dorsal >= Dv > dV >= ventral >= V)  \n  D:%i \n DV1:%i \n DV2:%i \n V:%i \n',DV(1),DV(2),DV(3),DV(4)),...
                        'Position',[0 yT+450 140 150],...
                        'HorizontalAlignment','Left');
                    imgSec = imgBase(:,:,zIx); 
                    hImg = image(imgSec);  axis equal; axis tight;
                    
                    if DV(3) > zIx && zIx > DV(2) || (zIx < DV(1) || DV(4)< zIx)  % this region will be deleted
                        set(h,'Color',[1 0 0]); 
                    elseif zIx == DV(1) || zIx == DV(2) || zIx ==DV(3) || zIx == DV(4) % boundary layer (included)
                        set(h,'Color',[0 0 1]); 
                    else
                        set(h,'Color',[1 1 1]*0.8); 
                    end
                    
                    k = 0;
                    btn = 0;
                    state = 'select spots';
                    while btn == 0
                        btn = waitforbuttonpress;
                        k = get(h,'CurrentCharacter');
                    end
                    switch lower(k)
                        case 'd'
                            zIx = zIx - 1;
                        case 'a'
                            zIx = zIx + 1;
                        case 'c'
                            zIx = zIx - 10;
                        case 'z'
                            zIx = zIx + 10;
                        case 's' % select spots
                            isSelSpots = 1;
                        case 'w'  % check again
                            xyz = floor([xyS'; zS]);                            
                            phiTheta = findThetaPhi(imgBase,xyz);
                            imgBase = correctLateralTilt(imgBase,phiTheta);  
                            zIx = DV(3);   
                            wSt = 1;        
                            description = 'tilt is corrected select DV regions';               
                            break
                        case 'q' % quit tilt correction
                            q = 1;
                            break
                        case 'e' % load the whole stack
                            if ~exist('lateralTilt.mat')
                                phiTheta = [0,0];
                                imgBase = img3D;
                                description = 'select spots to correct for the tilt';
                            else
                                load('lateralTilt.mat')
                                imgBase = correctLateralTilt(img3D,phiTheta); 
                                description = 'tilt is corrected select DV regions';
                            end
                            [~, mxind] = max(imgBase(:));
                            [~, ~, zIx] = ind2sub(size(imgBase),mxind);
                            zIx = DV(3);
                            wSt = 1;
                            break
                        case '1'
                            updDV = 1;
                            DV(1) = zIx; % dorsal
                            state = sprintf('dorsal crop selected: %d',zIx);
                        case '2' % update Dv value
                            updDV = 1;
                            DV(2) = zIx;
                            state = sprintf('DORSAL-ventral selected:%i',zIx);
                        case '3'
                            updDV = 1;
                            DV(3) = zIx;
                            state = sprintf('dorsal-VENTRAL selected:%i',zIx);
                        case '4' % update vD value
                            updDV = 1;
                            DV(4) = zIx; % ventral
                            state = sprintf('ventral crop selected: %d',zIx);
                    end
                    if nS>1
                    if lower(k) == 's' && (xyS(nS-1) == 0 || zS(nS-1) == 0)
                        zIx = zIx - 1;
                        nS = nS - 1;
                        q = 0;
                        display('wrong input detected, re-enter');
                    end
                    end
                    %% keep the z-scan in the valid range
                    if wSt == 1 %whole stack
                        if zIx > size(imgBase,3)
                            zIx = size(imgBase,3);
                        elseif zIx < 1
                            zIx = 1;
                        end
                    else
                        if zIx > fr2
                            zIx = fr2;
                        elseif zIx < fr1
                            zIx = fr1;
                        end
                    end
                end % while (1)
            end % while (quit)
            if updDV % DV layers updated
                display('DV updated')
                save('DV.mat','DV'); % =[D Dv dV V]
            end 
            if ~exist('phiTheta')
                phiTheta=[0 0];
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
    %% folders
    % 3D
    if 7~=exist('Rotated3D')
        mkdir('Rotated3D');
    end
    if 7~=exist('Rotated3D/3D')
        mkdir('Rotated3D/3D');
    end
    % dorsal-ventral 3D
    if 7~=exist('Rotated3D/dorsal3D')
        mkdir('Rotated3D/dorsal3D');
    end
    if 7~=exist('Rotated3D/ventral3D')
        mkdir('Rotated3D/ventral3D');
    end    
    %% cropping
    img3Dventral = img3D(:,:,DV(1):DV(2));
    img3Ddorsal = img3D(:,:,DV(3):DV(4));
    img3D = img3D(:,:,DV(1):DV(4));
    %% writing 3D stacks
    fnameRtd = strcat('Rotated3D/3D/',fNameWr);
    stackWrite(img3D,fnameRtd);
    fnameDorsal3D = strcat('Rotated3D/dorsal3D/',fNameWr);
    stackWrite(img3Ddorsal,fnameDorsal3D);
    fnameVentral3D = strcat('Rotated3D/ventral3D/',fNameWr);
    stackWrite(img3Dventral,fnameVentral3D);
    %% max-projection
    maxP        = max(img3D,[],3);
    maxPdorsal  = max(img3Ddorsal,[],3);
    maxPventral = max(img3Dventral,[],3);
    
    fnameMxProj         = strcat('Rotated3D/3D/','maxProj-',fNameWr);
    fnameMxProjDorsal   = strcat('Rotated3D/dorsal3D/','maxProj-',fNameWr);
    fnameMxProjVentral  = strcat('Rotated3D/ventral3D/','maxProj-',fNameWr);
    imwrite(uint16(maxP),char(fnameMxProj),'tiff');
    imwrite(uint16(maxPdorsal),char(fnameMxProjDorsal),'tiff');
    imwrite(uint16(maxPventral),char(fnameMxProjVentral),'tiff');
    

    
end

%[maxPventral MxZeInd]= max(img3D(:,:,DV:DV+5),[],3);  figure(1); imagesc(maxPventral); figure(2); imagesc(MxZeInd);
%[maxPventral MxZeInd]= max(img3D(1:100,1:100,DV:DV+5),[],3); figure(1); imagesc(maxPventral); figure(2); imagesc(MxZeInd);
