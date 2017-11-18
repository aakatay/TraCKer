%% display multiple images
if noTile || isDebugFilt 
    nFrTile = 1; % number of frames to be tiled
    nRow = 1;  
else
    mag = 10; % magnification
    [mag_, ~, ~, ~]= calcMaxMag(zeros(Boy1,En1),mag);
    if mag_<mag, mag=mag_; end;
    temp=get(0,'ScreenSize');
    screenX = temp(3); screenY = temp(4);
    nX = round((screenX + 200)/En1/mag);
    nY = round((screenY - 0)/Boy1/mag);

    nFrTile = nX*nY; % number of frames to be tiled
    nRow = nY;  
    if nFrTile > Frames, nFrTile =1;nRow=1;    end
end
    
% frm = 1;
framesVec =  1 +(nTile-1)*nFrTile : nTile*nFrTile;
for j=framesVec
    k = j - frstFrm -(nTile-1)*nFrTile + 1; % frame number in the tile
    for iAv = 1 : binFrame % bin frames
        frmRead = (j-1)*binFrame+iAv;
        temp1(:,:,iAv) = double(imread(fname,frmRead+frstFrm2-1));
        if isBALM, 
            frmRead = frmRead+1; 
            temp0(:,:,iAv) = double(imread(fname,frmRead+frstFrm2-2));
        end;
    end
    CoeffFitTile_(:,:,k) = ones(size(temp1))*CoeffFit(j);
    IMG = mean(temp1(yy1:yy2,xx1:xx2,:),3);
    if isBALM
        IMG0 = mean(temp0(yy1:yy2,xx1:xx2,:),3);
       IMG = IMG-IMG0;
    end
    IMGstack(:,:,k) = uint16(IMG);
    imgInRound = IMG; % no smoothing
    
    % DETECTION FILTER
    detectFilter; % IMG --> IMGfilt
    IMGfiltStack(:,:,k) = IMGfilt;
    
end

%nFr = size(IMGstack,3); % number of frames
nCol =  nFrTile/nRow;
nTileXY= [nRow,nCol];

% tile
CoeffFitTiled = tileFrames(CoeffFitTile_,nRow);
IMGtiled = tileFrames(IMGstack ,nRow);
IMGfiltTiled = tileFrames(IMGfiltStack ,nRow);
SHOW = IMGtiled;


