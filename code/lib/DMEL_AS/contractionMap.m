
    clear;
    close all;
    F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)

    Npx=128;
    Npx=512;
    Npx=201;
    minDist = 40; % min distance to evaluate trace speed
    imgFnm = 'maxProj.tif';
    if ~exist(imgFnm), imgFnm2 = rdir(sprintf('..\\%s',imgFnm)); imgFnm=imgFnm2.name; end;
    [yNpx,xNpx]=size(imread(imgFnm));

    dT = 1; % time interval [sec]
    BinSz = 8;
    %% load trace data
    d = rdir('traceData-coeff*.mat');
    [Y,J] = sort([d.datenum]); 
    traceDataPath = {d(J(end)).name}; % load the most recent posData-coeff
    load(char(traceDataPath));

    TraceX(isnan(TraceX)) = 0;
    TraceY(isnan(TraceY)) = 0;
    TraceINT(isnan(TraceINT)) = 0;
    X = TraceX;
    Y = TraceY;
    INT = TraceINT;
    INT(find(INT<0))=min(INT(find(INT>0)));
    
    dataCrop = 0;
    if dataCrop
        TraceX = TraceX - min(TraceX(:));
        TraceY = TraceY - min(TraceY(:));
    end
    
    %% colormap
    %colormap('jet'); 
    CM = colormap;close;
    intMax = max(INT(:));
    f = 1; % frame #
    R1 = 1; R2 = 10;
    C = CM(round(INT(R1:R2,f)/intMax/2*64+32),:);
    [Ntr T] = size(TraceX);

    fnameOut = '_roStack.tif';
    %figure; maximize
    %h = waitbar(0,'calculating the density map...');
    NbinsX = xNpx/BinSz;
    NbinsY = yNpx/BinSz;
    for t = 1:T-1 % # of frames  
        for i = 1:NbinsX % X
            for j = 1:NbinsY % Y
                findTracesApproach; % diffDr = findTracesApproach(TraceX,TraceY,T,i,j,BinSz)
                if mean(dRs)  < minDist
                    contrMap(i,j,t) = mean(approachParam)*40;
                else
                    %wait(0.1);
                end
                disp(t/T);
            end
        end
        subplot(1,2,1)
        imagesc(contrMap(:,:,t)); axis image;
        title(sprintf('frame#: %i',t));
        %waitforbuttonpress
        fname =  sprintf('%s_%02i.tif',fnameOut,t);
        imwrite(uint16(contrMap(:,:,t)),fname);

    end
    
    %% save data
    savedata = roStack;
    %savedata = pitStack;
    savedata = 2^8*savedata/max(savedata(:));
    for t = 1:size(savedata,3) % # of frames  
        if t == 1
            imwrite(uint16(savedata(:,:,t)),fnameOut);
        else
            imwrite(uint16(savedata(:,:,t)),fnameOut,'tiff','WriteMode','append');
        end
    end
    close(h)
    return

    [TrIx_,FrIx_]=ind2sub(size(TraceX),find(TraceX>0)); % trace and frame indices
