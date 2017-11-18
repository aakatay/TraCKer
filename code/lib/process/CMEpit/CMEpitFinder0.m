% using recruitment data find pits defined by size intensity and 0
% prebleach intensty
% INPUT:
% 1- acq_*** data file 
% 2- acqTL.tif time lapse movie (optional)
% 3- binImgRcrt** files

close all;
clear all;

isDyn = 0;

isFilterPreBleach = 1; % filter out the ones appear in prebleach image
isFilterPreBleach2 = 0; % filter out the ones recruits in first bin frame

isIllumMask = 1;

% pit selection
mxSz = 45;
mnSz = 5;
minInt = 10;

if isDyn % ow. AP2
    mnSz = 3;
    minInt = 5;
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
TLimgFN = 'TLimg.tif';
pitRecruitmentProfileFN = 'pitRecruitmentProfile.tif'; % rec time profiles
pitsResultsFN = 'pitResults.tif';

pitCoors = 'pitCoors.mat';
illumMaskMAT = 'illumMask.mat';
CoeffPreBleach = 'CoeffPreBleach.mat';
pitRecImg = 'pitRecImg.mat';
% color coding of struct type
% 1: no type
% 2: nonbleaching (hot spots)
% 3: low intensity
% 4: bona-fida pits
% 5: plaques
% 6: very small

%% prebleach image
pbImg = imread(acqFN,1); % pre bleach image
pbImg = pbImg(1:256,1:256);
imwrite(pbImg,preBleachFN);
if ~exist(illumFN), illumFN = preBleachFN; end

%%
c1 = [0 1 0]; % green
c2 = [1 1 0]; % yellow (non-bleaching pits)
c3 = [0 1 1]; % cyan (low intensity pits < 5recs)
c4 = [0 1 0]; % green (pits)
c5 = [1 0 1]; % magenta ( large structs >45)
c6 = [1 0 0]; % red (small structs <3px)
CMpitTypes = [0 0 0; c1; c2; c3; c4; c5; c6];

R = double(imread(binImgRcrtSumFN));
R = R(1:1024,1:1024);

x1 = 65;dx = 128;y1 = 65;dy = 128;

x1 = 257;dx = 512;y1 = 257;dy = 512;

x1 = 1;
dx = 1024;
y1 = 1;
dy = 1024;

x2 = x1+dx-1;
y2 = y1+dy-1;

R = R(y1:y2,x1:x2);


        % time lapse mean
        if exist(acqTLfn)
            timelapseInf = imfinfo(acqTLfn);
            for i = 1:numel(timelapseInf)
                acqTL = double(imread(acqTLfn,i));
                acqTL = acqTL(1:256,1:256);
            end
%            TLimg = mean(acqTL(129:384,129:384),3);
            TLimg = mean(acqTL,3);
            imwrite(TLimg,TLimgFN);
        else
            in=input('no time lapse movie available. OK ... ?','s');
            TLimg = zeros(256);
        end

        %% illumination area selection
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
        end
        if ~isIllumMask, illumMask = ones(size(illumMask)); end;
        %% structure detection
        cvSz=3;
        R = R.*illumMask;
        Rcv = conv2(R,ones(cvSz),'same');
        %Rcv(1,:)=9;Rcv(end,:)=9;Rcv(:,1)=9;Rcv(:,end)=9;
        RcvBin = im2bw(Rcv,0);
        RcvBininv = 1-RcvBin;
        


        % remove internal space
        [~,Lt00]=bwboundaries(Rcv,'noholes'); % tight boundaries
        Rcv2 = RcvBininv;
        Rcv2(Lt00>1)=0;
        Rcv3 = conv2(Rcv2,ones(cvSz),'same');
        Rcv3_ = Rcv3;
        Rcv3(Rcv3>0)=1;
        if rem(cvSz,2)==0, Rcv3=circshift(Rcv3,[1 1]); end
        Rcv4 = 1-Rcv3;

        [Bt0,Lt0]=bwboundaries(Rcv4,'noholes'); % tight boundaries
        if 0 
            figure(1001); imagesc(R); axis image;
            figure(1002); imagesc(Rcv); axis image;
            figure(1003); imagesc(RcvBin); axis image;
            figure(1005); imagesc(RcvBininv); axis image;
            figure(1006); imagesc(Rcv3_); axis image;
            figure(1007); imagesc(Rcv3); axis image;
            figure(1008); imagesc(Rcv4); axis image;
            figure(1010); plotBoundaries(Bt,sum(Lt,3).*R,1); axis image;
        end

        %% find the number of recruitments to pits
        Bt={};s=1;Lt=[];LtInt=[];
        Lt0col = uint16(im2bw(Lt0,0));
        ix0 = [];
        for i = 1:numel(Bt0) %find each tight section area
            typ=1;
            Lt0_ = zeros(size(Lt0));
            Lt0_(Lt0==i)=1;
            nr = sum(sum(Lt0_.*R)); % number of recruitments
            if nr <= 4, continue; end
            
            if sum(Lt0_(:)) > mxSz, typ=5; end
            if sum(Lt0_(:)) < mnSz, typ=6; end
            Lt0col(Lt0 ==i) = typ;
            if typ==5 || typ ==6, continue; end;

            Lt(:,:,s) = Lt0_;
            LtInt(:,:,s) = Lt0_*nr;
            Int(s) = nr;
            
            if nr <minInt    
                Lt0col(Lt0 ==i) = 3;
            end
            Bt = [Bt;Bt0(i)];
            s = s + 1;
            ix0 = [ix0 i];
        end
        ixPit = find(Int>=minInt);
        Lt = Lt(:,:,ixPit);
        Bt = Bt(ixPit);
        LtInt = LtInt(:,:,ixPit);
        Int = Int(ixPit);
        ix0 = ix0(ixPit);
        
        pitR = sum(Lt,3).*R;
        figure(1011); imagesc(pitR); axis image;
        
        pitSurf = sum(LtInt,3); 
        figure(1111); imagesc(pitSurf); axis image;
        colorbar;
        
        imwrite(uint16(pitR),recFN)
        imwrite(uint16(pitSurf),recIntFN)
        %% pick new started pits
        
        pbImg = double(imread(preBleachFN,1)); % pre bleach image
        
        %% find prebleach coeff
        % guess coeff from background
        if exist(CoeffPreBleach)
            load(CoeffPreBleach);
            [din] = CMEpitFinder_detectionLoop(pbImg,CoeffThresh);
            dinBW = im2bw(din,0);            
        else
            uiwait(msgbox('select background for prebleach image','Success','modal'));
            figure(121)
            imagesc(pbImg)
            [roi,r] = imcrop(gcf);
            bckgrnd = mean(roi(:));

            mx = max(pbImg(:));
            CoeffThresh=bckgrnd*3;
            [din] = CMEpitFinder_detectionLoop(pbImg,CoeffThresh);
            hFig=figure(100);
            pbImgFilt = pbImg - pbImg.*im2bw(din,0);

            pbImgFilt = pbImg;
            dinBW = im2bw(din,0);
            pbImgFilt(dinBW)=mx;
            imagesc(pbImgFilt)

            % coeff loop
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
                        break;
                end
                [din] = CMEpitFinder_detectionLoop(pbImg,CoeffThresh);
                delete(100)
                hFig=figure(100);
                pbImgFilt = pbImg;
                dinBW = im2bw(din,0);
                pbImgFilt(dinBW)=mx;
                imagesc(pbImgFilt)
            end
            save(CoeffPreBleach,'CoeffThresh');
        end
        imwrite(uint16(dinBW),prebleachStructsFN)
        %% remove pre bleach and fade structures   
        
        Rfrst5sec = double(imread(binImgRcrtFN,1));
        Rfrst5sec = Rfrst5sec(1:1024,1:1024);
        dinBW = repelem(dinBW,4,4);
        Lt2 = Lt;
        ixdel = [];
        for i =1:size(Lt2,3)
            po = Lt2(:,:,i)+isFilterPreBleach*im2bw(dinBW+isFilterPreBleach2*Rfrst5sec,0); % prebleach overlay
            if ~isempty(find(po>=2))
                ixdel=[ixdel i]; % remove overlapping structures
                Lt0col(Lt0 == ix0(i)) = 2; %nonbleaching
            else
                Lt0col(Lt0 == ix0(i)) = 4; % pits
            end
        end
        Lt(:,:,ixdel) = [];
        Int(ixdel) = [];
        LtInt(:,:,ixdel) = [];
        Bt(ixdel) = [];
        sz0 = size(Lt2,3);
        sz1 = size(Lt,3);
        perc = round((sz0-sz1)/sz0*100);
        disp(sprintf('%i%% of pits are in the prebleach image',perc))
        
        pitRfilt = sum(Lt,3).*R;
        pitSurf = sum(LtInt,3); clear LtInt;
        imwrite(uint16(pitRfilt),recPitsFN)
        imwrite(uint16(pitSurf),recPitsIntFN)
        
        Lt0colRGB = ind2rgb(uint16(Lt0col),CMpitTypes);
        imwrite(Lt0colRGB,recStructsFN);
        
        %% index pits
        A=zeros(size(R));
        figure(99)
        imagesc(A)
        pr=[]; % pit rec images
        szpad = 8;
        XY = [];
        for i =1:size(Lt,3)
            pr0 = Lt(:,:,i).*R;
            [Xc,Yc] = centOfMassLocCore(pr0,size(pr0,1));
            xy = [Xc,Yc];
            text(xy(1),xy(2),num2str(i))
            
            dxy = (sign(xy - round(xy)) + 1)/2;
            rxy = round(xy);
            rxy = rxy + dxy+szpad;
            x=uint16(rxy(1)); y=uint16(rxy(2));
            pr_ = padarray(pr0,[szpad szpad]);
            pr(:,:,i) = pr_(y-8:y+7,x-8:x+7);
            
            XY = [XY; round(xy)];
        end

        pos2 = [0 0 size(R,1) size(R,2)];
        set(gcf,'units','pixels','Position',pos2); 
        set(gca,'units','pixels','Position',pos2); 
        
        imgFig = getframe(gcf); 
        figCap_pitNumbers = imgFig.cdata;
        imwrite(figCap_pitNumbers,recPitNumbersFN) 
        save(pitCoors,'XY');
            
        
        %% find the duration of recruitment
        Rt = [];nr=[];
        iminf = imfinfo(binImgRcrtFN);
        nf = numel(iminf);
        
        for i = 1:nf
            Rt(:,:,i) = imread(binImgRcrtFN,i);
        end
        Rt = Rt(y1:y2,x1:x2,:);
        
        rp = zeros(nf,1); % recruitment profile (frames when recruited)
        strInt =[];
        for i =1:size(Lt,3)
            Rt_ = Rt.*repmat(Lt(:,:,i),1,1,nf);
            %if ismember(i,ixdel), continue; end;
            nr(:,i) = sum(reshape(Rt_,dy*dx,nf));
            nr(:,i) = nr(:,i)/max(nr(:,i))+i;
            strInt = sprintf('%s,''No%i,int:%i''',strInt,i,Int(i));
        end
        strInt(1) = [];
        
        
        %% print pitRecruitmentProfile
        figure(1000)
        maximize;
        t = (1:nf)*5;
        plot(t,nr,'Linewidth',2)
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

        pos2 = [0 0 size(R,1) size(R,2)];
        set(gcf,'units','pixels','Position',pos2); 
        set(gca,'units','pixels','Position',pos2); 
        imgFig = getframe(gcf); 
        figCap_recProf = imgFig.cdata;
        
        
        save(pitRecImg,'pr'); % pit rec images
        
        
        %% compile tif stack
        Rt = double(imread(binImgRcrtTimeFN));
        Rt = Rt(1:1024,1:1024);
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
                
        %3 acqPreBleach
        mx = max(pbImg(:));
        pbImg_ = repelem(double(pbImg)/mx*256,4 ,4);
        pbImg_  = ind2rgb( round(pbImg_) ,CMgray);
        imwrite(pbImg_,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        %4  pits_prebleachStructs
        dinBW_ = ind2rgb(dinBW*256,CMgray);
        imwrite(dinBW_,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        %5 pit types
        Lt0colRGB = ind2rgb(uint16(Lt0col),CMpitTypes);
        imwrite(Lt0colRGB,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        %6 pits
        mx = max(pitR(:));
        pitR_  = ind2rgb( uint16(pitR/mx*256) ,CMparula);
        imwrite(pitR_,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        %7 incoming pits
        mx = max(pitRfilt(:));
        pitRfilt_  = ind2rgb( uint16(pitRfilt/mx*256) ,CMparula);
        imwrite(pitRfilt_,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        %8 pit numbers
        figCap_pitNumbers_ = im2bw(figCap_pitNumbers);
        figCap_pitNumbers_  = ind2rgb( uint16(figCap_pitNumbers_*256) ,CMgray);
        imwrite(figCap_pitNumbers_,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        %9 pit time profiles
        imwrite(figCap_recProf,pitsResultsFN,'WriteMode','append','Compression', 'none') 
                         
        %10 time lapse mean
        mx = max(TLimg(:));
        TLimg_ = repelem(double(TLimg)/mx*256,4 ,4);
        TLimg_  = ind2rgb( round(TLimg_) ,CMgray);
        imwrite(TLimg_,pitsResultsFN,'WriteMode','append','Compression', 'none') 
        
        
        
        
        %1 recruitments time
        %2 recruitments
        %3 acqPreBleach
        %4 pits_prebleachStructs
        %5 pit types
        %6 recruitments for pits
        %7 recruitments for incoming pits
        %8 pit numbers
        %9 pit recruitment time profiles
        %10 time lapse mean (if there's one)