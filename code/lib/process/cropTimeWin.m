% crops recruimentTrack data in time
% run in recComp
close all;

cd ..
fimg = rdir('AVG_lap*'); % image file
binFN = 'binImgRcrt.tif';
binTimeFN = 'binImgRcrtTime.tif';
timeSelFN='recTimeSel.txt'; % for time selection
driftYesFN = 'driftDetected.txt';
driftNoFN = 'driftDetected_NO.txt';
avgFN = fimg.name;


if exist(timeSelFN)
    frt = fopen(timeSelFN);
    Rtxt = fgetl(frt);
    fclose(frt);
    recFrames =sscanf(Rtxt,'%i:%i'); % recBin frame no's
    f1 = recFrames(1);
    f2 = recFrames(2);
end

%% drift check
fdrift=rdir('driftDetected*.txt');
isCropTimeWin = 0;
if isempty(fdrift) % skip if done before
    isCropTimeWin = 1;
    fnameRecSumDIR = rdir('binImgRcrtSum_time*');
    fnameRecTimeDIR = rdir('binImgRcrtTime_time*');
    for i = 1:numel(fnameRecSumDIR)
        delete(fnameRecSumDIR(i).name)
        delete(fnameRecTimeDIR(i).name)
    end
    %% recTimeSel

    infBin = imfinfo(binFN);
    cropPos = [ 1 1 infBin(1).Width infBin(1).Height];
    nf = numel(infBin);

    if ~exist('f1')
        f1= 1;
        f2 = nf;
    end
    
    [f1,f2,Rt]=cropTimeWinCore(binFN,f1,f2,cropPos);
    
    frt = fopen(timeSelFN,'w');
    fprintf(frt,'%i:%i',f1,f2);
    clear binImg;

    %% drift correction
    sy = 0;
    sx = 0;
    while (1) 
        A0 = imread(avgFN);
        imagesc(A0);
        frt = fopen(timeSelFN);
        Rtxt = fgetl(frt);
        fclose(frt);

        recFrames =sscanf(Rtxt,'%i:%i'); % recBin frame no's
        recTimes = recFrames*5; % [seconds]
        ri = recTimes(1);
        rf = recTimes(2);
        ri = ri - 5;

        
        imagesc(A0);
        A = repelem(A0,4,4);
        % shift values for overlap correction
        A=circshift(A,sx,2);
        A=circshift(A,sy,1);

        Asc = double(A);
        Asc = ceil(Asc/max(Asc(:))*254);
        CM=gray(256);
        CM(end,:)=[1 0 0]; % red
        Asc(Rt>0) = 255;
        colormap(CM)
        imagesc(Asc)
        axis image
        maximize
        title(sprintf('sx:%i, sy:%i. use arrow keys to overlay coords. Press ''q'' to save&quit',sx,sy))
        btn = waitforbuttonpress;
        k = get(gcf,'CurrentCharacter');
        k = k + 0; % convert to number
        switch k 
            case 28 % left
                sx = sx - 1;
            case 29 % right
                sx = sx + 1;
            case 30 % up
                sy = sy - 1;
            case 31 % down
                sy = sy + 1;
            case 113 % q: quit
                break;
        end

    end
    imgOvFN = sprintf('recTimeSelOverlay_avgImage_%isx_%isy.tif',sx,sy);
    title(sprintf('Shifts needed for overlay. sx:%i, sy:%i',sx,sy))
    szx=size(Asc,2);    szy=size(Asc,1);
    set(gcf,'units','pixels','Position',[0,0,szx,szy+30]); 
    set(gca,'units','pixels','Position',[0,0,szx,szy]); axis tight;
    imgFig = getframe(gcf);
    dataImg = imgFig.cdata; 
    imwrite(uint16(dataImg),imgOvFN);

    if sx~=0 || sy~=0
        sprintf('Shifts needed for overlay. sx:%i, sy:%i',sx,sy)
        frt = fopen(driftYesFN,'w');
        fprintf(frt,'Shifts needed for overlay. sx:%i, sy:%i',sx,sy)
        fclose(frt)
        if exist(driftNoFN), delete(driftNoFN);end;
    else % no drift
        disp('No drift')
        frt = fopen(driftNoFN,'w');
        fprintf(frt,'Shifts needed for overlay. sx:%i, sy:%i',sx,sy);
        fclose(frt);

    end
    % write binImgRcrtSum with selected frames
    binImgSumFN = sprintf('binImgRcrtSum_time%i-%is.tif',ri,rf);
    binImgTimeFN = sprintf('binImgRcrtTime_time%i-%is.tif',ri,rf);
    binImg = [];
    j=1;
    for i = recFrames(1):recFrames(2)
        binImg(:,:,j) = imread(binFN,i);
        j = j +1;
    end
    
    %fid = fopen('frames.txt');
    %if fid>0, frstFrm = fscanf(fid,'%i:'); fclose(fid); else,frstFrm=1; end;
    
    %% time selected crops
    binImgSum = sum(binImg,3);
    imwrite(uint16(binImgSum),binImgSumFN);
    imwrite(uint16(Rt),binImgTimeFN);
    sxsy = [sx sy];
    figure;imagesc(im2bw(binImgSum,0));axis image
    figure;imagesc(im2bw(Rt,0));axis image
    
    ccc=21;
else
    if exist(driftYesFN)
        fdrift = fopen('driftDetected.txt');
        sxsy = fscanf(fdrift,'Shifts needed for overlay. sx:%i, sy:%i');
        imgOv = imread(sprintf('recTimeSelOverlay_avgImage_%isx_%isy.tif',sxsy(1),sxsy(2)));
        tittxt = sprintf('Shifts needed for overlay. sx:%i, sy:%i',sxsy(1),sxsy(2));
        fclose(fdrift);
    else
        sxsy = [0 0];
        imgOv = imread('recTimeSelOverlay_avgImage_0sx_0sy.tif');
        tittxt = 'no drift detected';
    end
    imshow(imgOv*2^8)
    title(tittxt);
    pause(0.2);
end
cd('recComp')




