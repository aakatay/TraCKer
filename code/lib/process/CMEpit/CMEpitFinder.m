% using recruitment data find pits defined by size intensity and 0
% prebleach intensty
% INPUT:
% 1- acq_*** data file 
% 2- acqTL.tif time lapse movie (optional)
% 3- binImgRcrt** files

close all;
clear all;

isDyn = 0;
isIllumMask = 0;
isSelectManual = 0;

% pit selection
mxSz = 45;
mxSz2 = 25;
minInt = 5;
minInt2 = 10;
smearThresh = 0.5;

if isDyn % ow. AP2
    minInt = 5;
    minInt2 = 10;
end

%% files
acqFN = rdir('acq_*'); % data file
acqFN = acqFN(1).name;

illumFN = 'acqIllum.tif';
binImgRcrtSumFN = 'binImgRcrtSum.tif';
binImgRcrtFN = 'binImgRcrt.tif';
binImgRcrtTimeFN = 'binImgRcrtTime.tif';
acqTLfn = 'acqTL.tif'; % time lapse file

% output files:
preBleachFN = 'acqPreBleach.tif'; % frst frame of acq*.tif
prebleachStructsFN = 'pits_prebleachStructs.tif';
recStructsFN = 'pits_recStructs.tif'; % color coding of struct type
recFN = 'pitRecruitments.tif';
recIntFN = 'pitIntensity.tif';
recPitsFN = 'pitRecruitmentsFresh.tif';
recPitsIntFN = 'pitIntensityFresh.tif';
recPitNumbersFN = 'pitNumbers.tif';
acqTLmeanFN = 'acqTLmean.tif';
pitRecruitmentProfileFN = 'pitRecruitmentProfile.tif'; % rec time profiles
pitsResultsFN = 'pitResults.tif';

pitCoors = 'pitCoors.mat';
illumMaskMAT = 'pit_IllumMask.mat';
CoeffPreBleach = 'pit_CoeffPreBleach.mat';
pitRecImg = 'pitRecImg.mat';
% color coding of struct type
% 1: no type
% 2: nonbleaching (hot spots)
% 3: low intensity
% 4: bona-fida pits
% 5: plaques
% 6: very small

%% prebleach image
fnFrstFrm = 'pitAcqFrstFrmNo.txt';
fid = fopen(fnFrstFrm);
if exist(fnFrstFrm), frstFrm = fscanf(fid,'%i',1); else, frstFrm = 1; end;
acqFrstFrm = double(imread(acqFN,frstFrm)); % pre bleach image
%if ~(size(acqFrstFrm,1)==szXY && size(acqFrstFrm,2)==szXY), error('size not 256by256'); end

szX = size(acqFrstFrm,2);
szY = size(acqFrstFrm,1);

imwrite(uint16(acqFrstFrm),preBleachFN);
if ~exist(illumFN), illumFN = preBleachFN; end
    

%% color coding of pit structures
c1 = [0 0 1];           % blue (low intensity)
c2 = [.61 .51 .74];     % purple (large structs)
% not selected
c3 = [0 1 1];           % cyan (non-bleaching pits)         in acqFrstFrm: yes | in rFrstFrm: yes
c4 = [1 0 1];           % magenta (in coming)                  in acqFrstFrm: no  | in rFrstFrm: yes
c5 = [0.91 0.41 0.17];  % orange (hotspot _ non-selected)   in acqFrstFrm: yes | in rFrstFrm: no
c6 = [1 0 0];           % red (new _ non-selected)         in acqFrstFrm: no  | in rFrstFrm: no
% selected
c7 = [1 1 0];           % yellow (hotspot _ selected)  
c8 = [0 1 0];           % green (new _ selected)  

% rest 
c9 = [1 1 1];           % white (single dots)

CMpitTypes = [0 0 0; c1; c2; c3; c4; c5; c6; c7; c8; c9];

%% time lapse mean
if exist(acqTLfn)
    timelapseInf = imfinfo(acqTLfn);
    for i = 1:numel(timelapseInf)
        acqTL(:,:,i) = double(imread(acqTLfn,i));
    end
    if size(acqTL,1)==512 & size(acqTL,2)==512
        acqTLmean = mean(acqTL(129:384,129:384),3); 
    elseif size(acqTL,1)==256 & size(acqTL,2)==256
        acqTLmean = mean(acqTL,3);
    elseif size(acqTL,1)==szY & size(acqTL,2)==szX
        acqTLmean = mean(acqTL,3);
    else % size doesnt match
        acqTLmean = zeros(szY,szX);
    end

    imwrite(acqTLmean,acqTLmeanFN);
else
    disp('no time lapse movie available. moving on ...');
    pause(3);
    acqTLmean = zeros(szY,szX);
end

%% illumination area selection
if isIllumMask
    if exist(illumMaskMAT)
        load(illumMaskMAT);
    else
        uiwait(msgbox('select illumination area','Success','modal'));
        ps = 100; % pad size
        A = imread(illumFN);
        imagesc(padarray(A,[ps ps])); axis image;
        elps= imellipse(gca);
        position = wait(elps);
        illumMask0 = createMask(elps); % mask
        illumMask = imcrop(illumMask0,[ps+1 ps+1 size(A,2)-1 size(A,1)-1]);
        illumMask = repelem(illumMask,4,4);
        save(illumMaskMAT,'illumMask');
        delete(1)
    end
else
    A = imread(illumFN);
    illumMask = ones(size(A)*4); 
end

%% prebleach image
% calculate coeff from background ROI
if exist(CoeffPreBleach)
    load(CoeffPreBleach);
    [din] = CMEpitFinder_detectionLoop(acqFrstFrm,CoeffThresh);
    dinBW = im2bw(din,0);       
    acqFrstFrmFiltBW = dinBW;     
else
    uiwait(msgbox('select background for prebleach image','Success','modal'));
    figure(121)
    imagesc(acqFrstFrm)
    [roi,r] = imcrop(gcf);
    bckgrnd = mean(roi(:));

    mx = max(acqFrstFrm(:));
    CoeffThresh=bckgrnd*3;
    [din] = CMEpitFinder_detectionLoop(acqFrstFrm,CoeffThresh);
    hFig=figure(100);
    acqFrstFrmFilt = acqFrstFrm - acqFrstFrm.*im2bw(din,0);

    acqFrstFrmFilt = acqFrstFrm;
    dinBW = im2bw(din,0);
    acqFrstFrmFilt(dinBW)=mx;
    imagesc(acqFrstFrmFilt)
    maximize;

    % coeff loop
    isBreak = 0;
    while 1
        text(-10,-10,sprintf('coeff:%i',round(CoeffThresh)))
        waitforbuttonpress
        k = get(hFig,'CurrentCharacter');
        switch lower(k)
            case 's'
                CoeffThresh =  CoeffThresh - 50; upd = 1;
            case 'd'
                CoeffThresh =  CoeffThresh + 50; upd = 1;
            case 'w'
                CoeffThresh =  CoeffThresh - 200; upd = 1;
            case 'e'
                CoeffThresh =  CoeffThresh + 200; upd = 1;
            case 'q'
                q = 1;
                figure(121);delete(121)
                figure(100);delete(100)
                isBreak = 1;
        end
        [din] = CMEpitFinder_detectionLoop(acqFrstFrm,CoeffThresh);
        figure(100);delete(100)
        hFig=figure(100);
        acqFrstFrmFilt = acqFrstFrm;
        dinBW = im2bw(din,0);
        acqFrstFrmFiltBW = dinBW;
        acqFrstFrmFilt(dinBW)=mx;
        imagesc(acqFrstFrmFilt)
        if isBreak, break; end;
    end
    save(CoeffPreBleach,'CoeffThresh');
end
imwrite(uint16(acqFrstFrmFiltBW),prebleachStructsFN)
imwrite(uint16(acqFrstFrmFiltBW),prebleachStructsFN,'WriteMode','append','Compression', 'none')
acqFrstFrmFiltBW = repelem(acqFrstFrmFiltBW,4,4);
        
if ~exist(binImgRcrtSumFN), disp('prebleach files are set'); return; end;

        %% recruitment data
        R = double(imread(binImgRcrtSumFN)); % recs
        Rt = double(imread(binImgRcrtTimeFN)); % time recs
        RR = [];
        iminf = imfinfo(binImgRcrtFN);
        nf = numel(iminf);
        for i = 1:nf
            RR(:,:,i) = imread(binImgRcrtFN,i);
        end

        %% recruitment image
        rFrstFrm = double(imread(binImgRcrtFN,1)); % recruitments in the first frame

        %% structure detection
        R = R.*illumMask;
        cvSz=3;
        [Bt0,Lt0] = findBoundariesCore(R,cvSz);

        %% calculate specs of the candidate pit 
        % OUT: (1)Lt: pit recruitments  (2)Lt0col: color coded pit types
        Bt={};s=1;Lt=[];LtInt=zeros(szY*4,szX*4);
        Lt0col = uint16(im2bw(Lt0,0));
            Lt = Lt0;
        ix0 = [];
        for i = 1:numel(Bt0) %find each tight section area
            Lt0_ = zeros(size(Lt0));
            Lt0_(Lt0==i)=1;
            
            % total number of recs
            nr_ = sum(sum(Lt0_.*R)); % number of recruitments
            if nr_ == 1, continue; end; % skip single points
            nr(s) = nr_;
            
            %size
            sz(s) = sum(Lt0_(:));
            
            % circularity
            stats = regionprops('struct',Lt0_,'MajorAxisLength','MinorAxisLength');
            length(s) = stats.MajorAxisLength; width(s)=stats.MinorAxisLength;
            
            %Lt(:,:,s) = Lt0_;
            %LtInt(:,:,s) = Lt0_*nr(s)/sz(s); % intensity image
            LtInt = LtInt + Lt0_*nr(s)/sz(s); % intensity image
            
            Bt = [Bt;Bt0(i)];
            s = s + 1;
            ix0 = [ix0 i]; % Lt mapping index
        end
        
        %% classify pits
        % low intensity
        ix = find(nr<minInt);
        Lt0col = im2bw(Lt0,0)*9;
        Lt0col(ismember(Lt0,ix0(ix))) = 1; 
        nonPit = ix;
        
        % large structs
        ix = find(sz>mxSz);
        Lt0col(ismember(Lt0,ix0(ix))) = 2; 
        nonPit = [ix nonPit];
        
        ixPit = [1:s-1]';
        ixPit(nonPit) = [];
                
        % PITS: 3:nonbleach, 4:incoming, 5:hotspot, 6:new
        for i =1:numel(ixPit)
            Lt0_ = zeros(size(Lt0));
            ixLt = ix0(ixPit(i));
            Lt0_(Lt0==ixLt)=1;
            pbOv = Lt0_+acqFrstFrmFiltBW; %  pre-bleach overlay
            rfOv = Lt0_+im2bw(rFrstFrm,0); %  recruitment first frame overlay
            isPbOv = 0;
            isRfOv = 0;
            if ~isempty(find(pbOv>=2)), isPbOv=1; end;
            if ~isempty(find(rfOv>=2)), isRfOv=1; end;
            
            %isRfOv = 0;
            if      isPbOv == 1 && isRfOv == 1 
                ixPit(i,2) = 3;
                Lt0col(Lt0 == ixLt) = 3;
            elseif  isPbOv == 0 && isRfOv == 1 
                ixPit(i,2) = 4;
                Lt0col(Lt0 == ixLt) = 4;
            elseif  isPbOv == 1 && isRfOv == 0 
                ixPit(i,2) = 5;
                Lt0col(Lt0 == ixLt) = 5;
            elseif  isPbOv == 0 && isRfOv == 0
                ixPit(i,2) = 6;
                Lt0col(Lt0 == ixLt) = 6;
            end
        end
        
        
        % prebleach ratio
        nPrebleachPits = sum(ixPit(:,2)==5) + sum(ixPit(:,2)==3);
        nPits = size(ixPit,1);
        perc = round(nPrebleachPits/nPits*100);
        disp(sprintf('%i%% of pits are in the prebleach image',perc))
        
        %% select pits
        % new pits
        ixSel_ = find(ixPit(:,2)==6); 
        ixSel = ixPit(ixSel_,1);
        
        % filter1: size and intensity
        ixSel2 = ixSel( find(nr(ixSel)>=minInt2 & sz(ixSel)<=mxSz2) ); 
        
        % filter2: time smear
        ixSel3 = [];
        for i = 1:numel(ixSel2)
            ix = ixSel2(i);
            ixLt = ix0(ix);
            Lt0_ = zeros(size(Lt0));
            Lt0_(Lt0==ixLt)=1;
            % smear
            rt = Rt.*Lt0_;            
            r = im2bw(R.*Lt0_,0);
            szXY = max(szX,szY);
            r2 = zeros(szXY*4);
            rt2 = r2;
            r2(1:szY*4,1:szX*4) = r;
            rt2(1:szY*4,1:szX*4) = rt;
            [Xc,Yc] = centOfMassLocCore(r2,szXY*4);
            [XcT,YcT] = centOfMassLocCore(rt2,szXY*4);
            ds = sqrt((XcT - Xc)^2+(YcT - Yc)^2); % shift due to smear
            L = length(ix);
            smearVal(i) = ds;
            if smearVal(i) < smearThresh
                ixSel3 = [ixSel3 ix];
                Lt0col(Lt0 == ixLt) = 8;
            end
        end
        Lt0colRGB = ind2rgb(uint16(Lt0col),CMpitTypes);
        
        
        RT = zeros(size(Lt0));
        for i = 1:numel(ixSel2)
            % scale color
            Lt0_ = zeros(size(Lt0));
            ixLt = ix0(ixSel2(i));
            Lt0_(Lt0==ixLt)=1;
            rt = Rt.*Lt0_;  
            mn = min(nonzeros(rt));
            rt = rt-mn+1;
            rt(rt<0)=0;
            %mx = max(rt(:));
            %rt = rt*255/mx;
            RT = RT + rt;
        end
        pitR0t = RT;
        figure;
        CM = colormap('parula');
        CM(1,:)=0;
        colormap(CM)
        imagesc(RT)
        axis image
        %% isSelectManual
        if exist('manualSelectROI.mat')
            load('manualSelectROI.mat');
        else
            ROI = [];
            if isSelectManual
                uiwait(msgbox('press ''s'' & select structures & press ''q''','Success','modal'));
                hFig = figure;
                imagesc(Lt0colRGB); maximize; axis image;
                while(1)
                    waitforbuttonpress
                    k = get(hFig,'CurrentCharacter');
                    switch lower(k)
                        case 's' % select
                            [roi,r] = imcrop(gcf);
                            ROI = [ROI; r];
                        case 'q' % quit
                            break;
                    end
                end
            end
            save('manualSelectROI.mat','ROI')
        end
        structs = Lt0col;
        structs(structs<3) = 0;
        structs(structs==9) = 0;
        structs = im2bw(structs,0).*Lt;
        structIx = [];
        nrManual = [];
        for i = 1:size(ROI,1)
            roi = round(ROI(i,:));
            x = [roi(1) roi(1) roi(1)+roi(3) roi(1)+roi(3) roi(1)];
            y = [roi(2) roi(2)+roi(4) roi(2)+roi(4) roi(2) roi(2)];
            roi = poly2mask(x,y,szY*4,szX*4);
            %imagesc(roi)
            stIx = unique(roi.*structs);
            stIx(1) = [];
            if numel(stIx)>1, warning('two structs selected'); break;end;
            structIx = [structIx stIx];
            nrManual(i) = sum(sum(roi.*R));
        end
        
        %% index pits & find the recruitment profile
        A=zeros(size(R));
        figure(99)
        imagesc(A)
        pr=[]; % pit rec images
        szpad = 8;
        XY = [];
        rp=[];
        strInt =[];
        LtPits = zeros(size(Lt));
        ixLT = [ix0(ixSel3) structIx];
        
        nrsel = [nr(ixSel3) nrManual];
        for i =1:numel(ixLT)
        %for i =1:numel(ixSel3)
            %ix = ixSel3(i);
            %ixLt = ix0(ix);
            ixLt = ixLT(i);
            Lt0_ = zeros(size(Lt0));
            Lt0_(Lt0==ixLt)=1;
            
            % index pits
            pr0 = Lt0_.*R;
            pr02 = zeros(szXY*4);
            pr02(1:szY*4,1:szX*4) = pr0;
            [Xc,Yc] = centOfMassLocCore(pr02,szXY*4);
            xy = [Xc,Yc];
            text(xy(1),xy(2),num2str(i))
            
            dxy = (sign(xy - round(xy)) + 1)/2;
            rxy = round(xy);
            rxy = rxy + dxy+szpad;
            x=uint16(rxy(1)); y=uint16(rxy(2));
            pr_ = padarray(pr0,[szpad szpad]);
            pr(:,:,i) = pr_(y-8:y+7,x-8:x+7);
            
            XY = [XY; xy];
            
            % find the recruitment profile
            Rt_ = RR.*repmat(Lt0_,1,1,nf);
            rp(:,i) = sum(reshape(Rt_,szX*szY*16,nf));
            rp(:,i) = rp(:,i)/max(rp(:,i))+i;
            strInt = sprintf('%s,''No%i,int:%i''',strInt,i,nrsel(i));
            
            % pit surface
            LtPits_ = zeros(size(Lt));
            LtPits_(Lt == ix0(ix))= 1;
            LtPits(:,:,i) = LtPits_;
        end
        pos2 = [0 0 size(R,2) size(R,1)];
        set(gcf,'units','pixels','Position',pos2); 
        set(gca,'units','pixels','Position',pos2); 
        
        imgFig = getframe(gcf); 
        figCap_pitNumbers = imgFig.cdata;
        imwrite(figCap_pitNumbers,recPitNumbersFN) 
        save(pitCoors,'XY','LtPits');
        
        strInt(1) = [];

        %% display and write
        rOther = zeros(size(Lt0col));
        pitRfilt = zeros(size(Lt0col));
        pitR = zeros(size(Lt0col));
        
        rOther(Lt0col==9) = 1;
        pitRfilt(ismember(Lt0col,[3:6])) = 1;
        pitR(ismember(Lt0col,[7:8])) = 1;
        
        rOther = rOther.*R;
        pitRfilt = pitRfilt.*R;
        pitR = pitR.*R;
        
        figure(1011); imagesc(pitR); axis image;
        LtInt = sum(LtInt,3); 
        figure(1111); imagesc(LtInt); axis image;
        colorbar;
        
        
        imwrite(uint16(pitRfilt),recPitsFN)
        imwrite(uint16(LtInt),recPitsIntFN)
        imwrite(Lt0colRGB,recStructsFN);
        imwrite(uint16(pitR),recFN)
        imwrite(uint16(LtInt),recIntFN)
        
        %% print pitRecruitmentProfile
        figure(1000)
        maximize;
        t = (1:nf)*5;
        plot(t,rp,'Linewidth',2)
        legendCall = sprintf('legend(%s)',strInt);
        eval(legendCall)

        title({'recruitment number vs time';'legend names are total number of recruitments for each pit'})
        xlabel('time (sec)')
        ylabel('# recruitments')
        
        posScreen2 = [1921        -423        1080        1853];
        set(gcf,'units','pixels','Position',posScreen2); 
        posScreen2_ = [50        50        1000        1753];
        set(gca,'units','pixels','Position',posScreen2_); 
        imgFig = getframe(gcf); 
        figCap_recProf = imgFig.cdata;
        imwrite(figCap_recProf,pitRecruitmentProfileFN,'Compression', 'none') 

        pos2 = [0 0 size(R,2) size(R,1)];
        set(gcf,'units','pixels','Position',pos2); 
        set(gca,'units','pixels','Position',pos2); 
        imgFig = getframe(gcf); 
        figCap_recProf = imgFig.cdata;
        save(pitRecImg,'pr'); % pit rec images
        
        %% compile tif stack
        delete(pitsResultsFN)
        CMparula = parula(256);
        CMgray = gray(256);
        CMcool = cool(256);
        CMcool(1,:)=0;
        CMparula(1,:)=0;
        %1 recs time
        mx = max(Rt(:));
        Rt_  = ind2rgb( uint16(Rt/mx*256) ,CMcool);
        imwrite(Rt_,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        %2 recs
        mx = max(R(:));
        R_  = ind2rgb( uint16(R/mx*256) ,CMparula);
        imwrite(R_,pitsResultsFN,'WriteMode','append','Compression', 'none') 
                
        %3 acqFrstFrm
        mx = max(acqFrstFrm(:));
        acqFrstFrm__ = repelem(double(acqFrstFrm)/mx*256,4 ,4);
        acqFrstFrm_  = ind2rgb( round(acqFrstFrm__) ,CMgray);
        imwrite(acqFrstFrm_,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        %4  pits_prebleachStructs
        acqFrstFrmFiltBW_ = ind2rgb(acqFrstFrmFiltBW*256,CMgray);
        imwrite(acqFrstFrmFiltBW_,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        %5 pit types
        Lt0colRGB = ind2rgb(uint16(Lt0col),CMpitTypes);
        imwrite(Lt0colRGB,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        %6 other recs
        mx = max(rOther(:));
        rOther_  = ind2rgb( uint16(rOther/mx*256) ,CMparula);
        imwrite(rOther_,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        %7 not selected pits
        mx = max(pitRfilt(:));
        pitRfilt_  = ind2rgb( uint16(pitRfilt/mx*256) ,CMparula);
        imwrite(pitRfilt_,pitsResultsFN,'WriteMode','append','Compression', 'none') 

        %9 selected pits - recruitment intensity
        mx = max(pitR0t(:));
        pitR0t_  = ind2rgb( uint16(pitR0t/mx*256) ,CMparula);
        imwrite(pitR0t_,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        %9 selected pits - recruitment times
        mx = max(pitR(:));
        pitR__ = uint16(pitR/mx*256);
        pitR_  = ind2rgb( pitR__ ,CMparula);
        imwrite(pitR_,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        %10 acqFrstFrm overlay
        ffOverlay = acqFrstFrm__/256;
        ffOverlay(:,:,2) = pitR__;
        ffOverlay(:,:,3) = 0;
        imwrite(ffOverlay,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        %11 pit numbers
        figCap_pitNumbers_ = im2bw(figCap_pitNumbers);
        figCap_pitNumbers_  = ind2rgb( uint16(figCap_pitNumbers_*256) ,CMgray);
        imwrite(figCap_pitNumbers_,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        %12 pit time profiles
        imwrite(figCap_recProf,pitsResultsFN,'WriteMode','append','Compression', 'none') 
                         
        
        %13 time lapse mean
        mx = max(acqTLmean(:));
        acqTLmean_ = repelem(double(acqTLmean)/mx*256,4 ,4);
        acqTLmean_  = ind2rgb( round(acqTLmean_) ,CMgray);
        imwrite(acqTLmean_,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        
        
        
        %1 recruitment times
        %2 recruitments
        %3 acqFrstFrm
        %4 pits_prebleachStructs
        %5 structure types
        %6  single events
        %7 recruitments for structures  (not selected)
        %9 recruitments for pits (selected for intensity and size)
        %9 recruitment times for pits (selected for intensity, size)
        %10 acqFrstFrm overlay
        %11 pit numbers for time smear selected pits
        %12 pit recruitment time profiles
        %13 time lapse mean (if there's one)