function rtDetectThresh
%% CORRECTION: remove a line cropping an area
% detects threshold for tracking
% output Coeff.mat
% RUN in waSeq\


% image =======
% filtered image =======
% filtered image peak MASK 
% image PEAKs =======
% histogram peak (image PEAKs)

%cd('E:\MATLAB\TIRFcalibration\data\Ata01_5_125X100Y50x50_realtime');    
    cd waSeq
    %test = []; save test test
    
    dbg = 0;

    load ..\cfgRT;
    acqTime = cfg.acqTime; % [s]
    ndigit = cfg.ndigit; % # of digits for sequence number
    w = cfg.w;
    h = cfg.h;
    bc = cfg.bc; % [nM] background concentration
    
    tic; logFN = cfg.logThresh; fid = fopen(logFN,'w'); wait = 0;
    clck = clock; fprintf(fid,'start time m= %2i secs=%6.03f\n',clck(5),clck(6));
    

    digitFormat = sprintf('%%0%1ii',ndigit);
    outDIR = 'tracker\';
    
    CoeffFN = [outDIR 'Coeff_' cfg.label '.mat'];
    Coeff = []; save(CoeffFN,'Coeff');

    % coefficient multiplier
    cm1 = bc/10;
    cm2 = bc/27;
    coeffMul1 = 1+0.4; 
    coeffMul2 = 1-0.15; 
    coeffMul1 = 1+cm1; % corrects shift of histogram peak to the right
    coeffMul2 = 1-cm2; % corrects for increased intensity fluctuations

    coeffMul = coeffMul1*coeffMul2; % 4nM background


    %% set filenames and array
    load('..\fname0.mat'); % % filename_WA_
    
    %% feedback to calling function
    fcall = 'rtDetectThresh';
    fdbck.inWait = 0;
    fdbck.inWaitCounting = 0;
    fdbck.inPause = 0;
    fdbck.inSave = 0;
    fdbck.inSaveCounting = 0;
    fdbck.inSaveCountingIX = 0;
    fdbck.inSaveCountingMAX = cfg.inSaveCountingMAX;
    fdbck.inStop = 0;    

    n = 1;
    while (1)
        time = toc; fprintf(fid,'while loop n=%3i time=%6.03f\n',n,time);

        while (1) % wait for update
            fnameSeq = [fname0 'WA_' num2str(n,digitFormat) '.tif'];
            if ~exist(fnameSeq) % wait for update
                if wait == 0, time = toc; fprintf(fid,'wait for   n=%3i time=%6.03f\n',n,time); wait = 1; end
                %return;
                [fdbck] = funcFeedback(cfg.msgTXT,fdbck,fcall);
                if fdbck.inStop, break;  end % STOP
            else
                time = toc; wait = 0; fprintf(fid,'updated    n=%3i time=%6.03f\n',n,time);
                break; % continue
            end 
        end

        %% image
        IMG = imread(fnameSeq,1);
%                    IMG(90:94,1:4)=IMG(95:99,1:4); % ======remove
        [~,c] = loghist(IMG,30,0);
        hstep = c(2)-c(1); % hist step

        if dbg
            nc = 6; % # columns
            figure(505);
            subplot(2,nc,1)
            imagesc(invert(IMG));
            colormap('gray');axis image; title('image')

            subplot(2,6,nc+1)
            [~,c] = loghist(IMG,30,1);
        end

        %% filtered image 
        gausKernelSz = 3;
        gausKernelSg = 0.7;
        gausKernelSzPadded = 5;

        pdSz = (gausKernelSzPadded - gausKernelSz)/2; % pad size
        gausHat1=fspecial('gaussian', gausKernelSz, gausKernelSg);
        lap=[-1,-1,-1;-1,8,-1;-1,-1,-1];

        IMGfilt2 = nan(size(IMG));
        IMGfilt1 = imfilter(IMG,gausHat1,'symmetric');
        IMGfilt1 = imfilter(IMGfilt1,lap,'symmetric');
        IMGfilt = max([IMGfilt1(:) IMGfilt2(:)],[],2);
        IMGfilt = reshape(IMGfilt,size(IMGfilt1));
        if dbg
            subplot(2,nc,2)
            imagesc(invert(IMGfilt))
            colormap('gray');axis image;title('filtered image')
            imwrite(IMGfilt,'IMGfilt.tif')
            subplot(2,6,nc+2)
        end

        [~,c] = loghist(IMGfilt,30,dbg);

        %% filtered image peak MASK 
        BW = imregionalmax(IMGfilt, 8);

        %% image PEAKs
        IMGpeak = uint16(BW).*IMG;
        if dbg
            subplot(2,nc,3)
            imagesc(invert(IMGpeak))
            colormap('gray');axis image; title('image PEAKs')
            imwrite(IMGpeak,'IMGpeak.tif')

            subplot(2,6,nc+3)
            [hc,c] = loghist(IMGpeak,50,dbg); % hist count and centers
            xlabel('intensity'); ylabel('counts')
        end
        %[hc,c] = loghist(IMGpeak,50,1); % hist count and centers


        %%  histogram peak(image PEAKs)
        [hc,c] = loghist(IMGpeak,50,0); % hist count and centers
        [~,locs] = findpeaks(hc);
        peakInt1 = c(locs(1));

        Coeff(n) = peakInt1*coeffMul;
        save(CoeffFN, 'Coeff');
        n = n + 1; % # of frames
        
    end
    n = n - 1;
end