        load spotSel; % LOADS parameters
        
        % gaus fit PSF parameters (LOADED)
        %intP1G = 5; intP2G = 5; % needs to be odd (PSF power)
                
        %% gaus fit starting values and constraint param.
        % ==2G fit==
        % peak
        sg = 0.73; %[px] sg = 0.73nm (FWHM=183nm)
        sr = 1+0.2; % tolerance multiplier 
        %tP2G_1 = 1.75; % tolerance position (LOADED)
        
        % background
        sg0min = 3;
        sg0 = 5; % starting guess
        sg0max = inf;
        % tP2G_2 = inf; (LOADED)       
                
        % ==intensity 1G and 2G==
        %minmaxIntTol = inf; (LOADED) 
        tol = minmaxIntTol + 1; % peak
        tol0 = tol; % background 

        % define output XY array
        [iyTrack, ixTrack] = find(xBW>0); % positions of spot data in the X array
        Xgaus = zeros(size(xBW));
        Ygaus = zeros(size(xBW));

        % crop 
        x3x1=4;y3y1=x3x1;x3x2=6;y3y2=x3x2;
        x5x1=3;y5y1=x5x1;x5x2=7;y5y2=x5x2;
        x7x1=2;y7y1=x7x1;x7x2=8;y7y2=x7x2;
        Window7by7 = spotWin(y7y1:y7y2,x7x1:x7x2,1);

        %define fit parameters
        [n,m]=size(Window7by7); %assumes that I is a nxm matrix
        px0 = double(ceil(m/2)); % center position 
        m2 = m + 2; n2 = n + 2; 
                



        % gaus fit grid
        [XX,YY]=meshgrid(1:n,1:m);%your x-y coordinates
        x(:,1)=XX(:); % x= first column
        x(:,2)=YY(:); % y= second column
        
        % gaus display grid
        [n2,m2]=size(spotWin(:,:,1)); %assumes that I is a nxm matrix
        [XX,YY]=meshgrid(1:n2,1:m2); %your x-y coordinates
        xSHOW(:,1)=XX(:); % x= first column
        xSHOW(:,2)=YY(:); % y= second column
        
        h = waitbar(0,'Gaus fit localization ....');
        nSpots = size(spotWin,3);
        Frames = size(xBW,2);
        nSpots = 200;
        frstSpot =1;
        errCenter = cell(1,nSpots);
        nSpotsDisp = 1000;
        G0reg = zeros(31,nSpotsDisp);
        G0nonReg = zeros(41,nSpotsDisp);











for i = frstSpot:nSpots % every spot 
    if i>1000, i2=nSpotsDisp+1; else i2 = i; end;
    %% min & max values & normalization & background removal
    Window9by9 = double(spotWin(:,:,i));
    Window7by7 = Window9by9(y7y1:y7y2,x7x1:x7x2);

    % find max pixel pos.
    [pixMaxY_,pixMaxX_] = find(Window7by7==max(Window7by7(:)));
    R_ = (pixMaxY_-4).^2+(pixMaxX_-4).^2;
    [~, ix] = sort(R_);
    pixMaxY(i)=pixMaxY_(ix(1)); pixMaxX(i)=pixMaxX_(ix(1));

    % find background 
    gausKernelSz = 5;
    gausKernelSg = 2.2;
    gausBlur=fspecial('gaussian', gausKernelSz, gausKernelSg);
    win9by9blur_ = imfilter(Window9by9,gausBlur,'symmetric');
    data9by9 = Window9by9 - win9by9blur_; % remove the background
    IntFlatBckgrnd(i) = min(min(data9by9(4:6,4:6)));
    data9by9 = data9by9 - IntFlatBckgrnd(i); % no negative values allowed

    % background position
    [yBckgrnd_,xBckgrnd_] = find(win9by9blur_==max(win9by9blur_(:)));
    R_ = (yBckgrnd_-5).^2+(xBckgrnd_-5).^2;
    [vv, ix] = sort(R_);
    yBckgrnd(i)=yBckgrnd_(ix(1)); xBckgrnd(i)=xBckgrnd_(ix(1)); 
    % peak position
    intPeak(i) = max(max(data9by9(4:6,4:6))); % max in 3by 3 window
    data9by9blur(:,:,i2) = win9by9blur_/intPeak(i);
    data9by9 = data9by9/intPeak(i); % spot peak normalized to one
    data7by7norm = data9by9(2:8,2:8); % normalized
    % background & positive data
    data7by7 = data7by7norm - IntFlatBckgrnd(i); % IntFlatBckgrnd is always negative
    % figure(1); imagesc([data9by9+data9by9blur(:,:,i2) zeros(9,2) data9by9 zeros(9,2) data9by9blur(:,:,i2)]); axis image; continue            
    %% read spot intensities
    spReg = spReg7by7(spReg7by7(:,i)>0,i);
    spNonReg = spNonReg7by7(spNonReg7by7(:,i)>0,i);

    nSptsReg = numel(spReg);
    nSptsNonReg = numel(spNonReg);
    nSpts(:,i) = [nSptsReg nSptsNonReg];            

    %% read the registered spot contribution: G0reg & fit7by7reg
    fit9by9reg_ = zeros(m2);
    nSptsReg = numel(spReg);
    nonRegDispOffset = [0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0]';
    if ~isempty(spReg) % define G0reg for calc. the image of pre-regis. spots
        G0reg(1,i)=0; % flat background
        IntReg = INTnorm(spReg)/intPeak(i);
        xReg = Xgaus(spReg)+px0-double(xBW(i)); % use pre fit values
        yReg = Ygaus(spReg)+px0-double(yBW(i)); 
        sgx = fitVal(4,spReg);
        sgy = fitVal(6,spReg);
        % spots
        for j =  1:nSptsReg 
            if j>nSptsReg, IntReg(j)=0; xReg(j)=0; yReg(j)=0; end;
            G0reg([2:6]+(j-1)*5,i) = [IntReg(j) xReg(j) sgx(j) yReg(j) sgy(j)]; % fit values of regjstered spots
        end
        G0regDisp = double(G0reg(:,i))+nonRegDispOffset(1:size(G0reg,1));
        fit9by9reg_ = funSpReg(G0regDisp,xSHOW);
        fit9by9reg_ = reshape(fit9by9reg_,m2,m2);
    end
    fit9by9reg(:,:,i2) = fit9by9reg_;
    fit7by7reg = fit9by9reg(2:8,2:8,i2);

    %% calculate flatbackground contribution from registered spots
    IntFlatBckgrndDiff9by9_ = zeros(m2);
    IntFlatBckgrndDiff9by9sum = zeros(m2);
    IntFlatBckgrndDiff9by9sum_ = zeros(m2);
    for j = 1:numel(spReg)
        IntFlatBckgrndDiff9by9_ = zeros(m2);
        fit9by9reg_ = fit9by9reg(:,:,i2);
        j_ = spReg(j);

        IntFlatBckgrndDiff = ( IntFlatBckgrnd(i)*intPeak(i) - IntFlatBckgrnd(j_)*intPeak(j_) ) / intPeak(i); % background difference between two fits
        x_ = xBW(j_)+px0-double(xBW(i))+1;
        y_ = yBW(j_)+px0-double(yBW(i))+1;
        x1 = x_-1; x2 = x_+1;
        y1 = y_-1; y2 = y_+1;
%                 if x1 < 1, x1 = 1; end; if y1 < 1, y1 = 1; end;
%                 if x2 > m, x2 = m; end; if y2 > m, y2 = m; end;
        % select pixels 
        cropROI = zeros(m2);
        cropROI(y1:y2,x1:x2) = 1;
        fit9by9reg_ = cropROI.*fit9by9reg_;
         % intensity threshold to select pixels to add background difference
        intThr = max(fit9by9reg_(:))/5;
        pxSel2addBckDiff = find(fit9by9reg_>intThr);
        IntFlatBckgrndDiff9by9_(pxSel2addBckDiff) = IntFlatBckgrndDiff; 
        IntFlatBckgrndDiff9by9sum = IntFlatBckgrndDiff9by9sum + IntFlatBckgrndDiff9by9_;
        IntFlatBckgrndDiff9by9sum_ = IntFlatBckgrndDiff9by9sum_ + double(IntFlatBckgrndDiff9by9_>0);
    end
    IntFlatBckgrndDiff9by9sum_(IntFlatBckgrndDiff9by9sum_==0)=1;
    IntFlatBckgrndDiff9by9sum = IntFlatBckgrndDiff9by9sum./IntFlatBckgrndDiff9by9sum_;
    IntFlatBckgrndDiff9by9(:,:,i2) = IntFlatBckgrndDiff9by9sum;
    IntFlatBckgrndDiff7by7 = IntFlatBckgrndDiff9by9(2:8,2:8,i2);
    % input data to the fit function: data7by7diff
    data7by7diff = data7by7 - fit7by7reg - IntFlatBckgrndDiff7by7; % background and pre-fit spots are substracted
    %%
    clear CLim;
    figure(1)
    subplot(2,2,1)
    imagesc(data7by7);
    gcaImg(1) = gca; title('data7by7')
    subplot(2,2,2)
    imagesc(fit7by7reg);
    gcaImg(2) = gca; title('fit7by7reg')
    subplot(2,2,3)
    imagesc(data7by7diff);
    gcaImg(3) = gca; title('data7by7diff')
    subplot(2,2,4)
    imagesc(IntFlatBckgrndDiff7by7);
    gcaImg(4) = gca; title(sprintf('IntFlatBckgrndDiff7by7:%.02f',IntFlatBckgrndDiff7by7(1)))
    colorbar
    for ii = [1:4]
        CLim(:,ii) = get(gcaImg(ii),'CLim');
    end
    CLimMax = max(CLim(2,:));
    CLimMin = min(CLim(1,:));
    CLim = [CLimMin CLimMax];
    for ii = [1:4]
        set(gcaImg(ii),'CLim',CLim);
    end            
    %%

    %% define fit function and constraints for all possible spots
    % 8 spots :  1: bckgrnd 2: center spot 3-8: neigh. spots
    ySpt = px0;
    xSpt = px0;
    IntNonReg = zeros(6,1); % neighbouring spots
    xNonReg = zeros(6,1);
    yNonReg = zeros(6,1);

    if nSptsNonReg>0
        xNonReg(1:nSptsNonReg) = xBW(spNonReg)-xBW(i)+px0;
        yNonReg(1:nSptsNonReg) = yBW(spNonReg)-yBW(i)+px0;
    end
    IntBckgrnd = data9by9blur(yBckgrnd(i),xBckgrnd(i),i2);
    % position tolerances
    tP2G_1 = 1; % main peak(s)
    tP2G_2 = 2; % secondary peak (background)
%             G0nonReg(7:11,i2) = [ IntBckgrnd xBckgrnd(i) sg0 yBckgrnd(i) sg0]; % background 
%             LB(7:11,i2) = [0 0 xBckgrnd(i)-tP2G_2 sg0min yBckgrnd(i)-tP2G_2 sg0min];
%             UB(7:11,i2) = [ IntBckgrnd*tol0 xBckgrnd(i)+tP2G_2 sg0max yBckgrnd(i)+tP2G_2 sg0max];
    G0nonReg(1:6,i2) = [IntFlatBckgrnd(i) 1 xSpt sg ySpt sg]; % main spot
    LB(1:6,i2) = [0 1/tol xSpt-tP2G_1 sg/sr ySpt-tP2G_1 sg/sr];
    UB(1:6,i2) = [inf 1*tol xSpt+tP2G_1 sg*sr ySpt+tP2G_1 sg*sr];
    LB=double(LB);UB=double(UB);            
    % neighbouring spots
    %spNonRegMax = 6; % 6 spots (neighbouring spots)
    for j =  1:nSptsNonReg
        j_ = spNonReg(j);
        x1 = xBW(j_)-xBW(i)+px0-1; x2 = xBW(j_)-xBW(i)+px0+1; 
        y1 = yBW(j_)-yBW(i)+px0-1; y2 = yBW(j_)-yBW(i)+px0+1; 
        if x1 < 1, x1 = 1; end; if y1 < 1, y1 = 1; end;
        if x2 > m, x2 = m; end; if y2 > m, y2 = m; end;                
        crop3by3NonReg = data7by7norm(y1:y2,x1:x2);
        IntNonReg(j) = max(max(crop3by3NonReg)); 
        G0nonReg([12:16]-5+(j-1)*5,i2) = [IntNonReg(j) xNonReg(j) sg yNonReg(j) sg]; % fjrst guess values for the fjt
        LB([12:16]-5+(j-1)*5,i2) = [IntNonReg(j)/tol xNonReg(j)-tP2G_1 sg/sr yNonReg(j)-tP2G_1 sg/sr];
        UB([12:16]-5+(j-1)*5,i2) = [IntNonReg(j)*tol xNonReg(j)+tP2G_1 sg*sr yNonReg(j)+tP2G_1 sg*sr];            
    end

    %% fit function: funSpNonReg
    if ~exist('txFitSpNonRegFunc')
        [funText, txFitSpNonRegFunc] = genFitFunc(nSptsNonReg); % function text
    else
        [funText, ~] = genFitFunc(nSptsNonReg,txFitSpNonRegFunc); % function text
    end
    eval(funText);

    %% 3by3 and 2by2 data
    pxYc = 4; pxXc = 4;
    data3by3 = data7by7(pxYc-1:pxYc+1,pxXc-1:pxXc+1); 
    % find the 4 peak pixels
    dat=data3by3;
    clear ixX ixY;
    datCenter = dat(2,2); dat(2,2)=0;
    [~, iX]=max(dat(:)); dat(iX) = 0;
    ixY(1) = 2; ixX(1) = 2;
    [ixY(2), ixX(2)]=ind2sub(size(dat),iX);
    [~, iX]=max(dat(:)); dat(iX) = 0;
    [ixY(3), ixX(3)]=ind2sub(size(dat),iX);

    % find coor of 2by2 peak
    iserr2by2_ = [];
    if round(mean(ixX)) == max(ixX), ixX(4) = min(ixX); elseif round(mean(ixX)) == min(ixX), ixX(4) = max(ixX); else, iserr2by2_='split peak'; end
    if round(mean(ixY)) == max(ixY), ixY(4) = min(ixY); elseif round(mean(ixY)) == min(ixY), ixY(4) = max(ixY); else, iserr2by2_='split peak'; end

    if ~isempty(iserr2by2_) % split peaks
        dat=data3by3;
        if ixY(2)== 2 && ixX(2) == 1 % cases 2 & 4
            ixX(3:4) = [1,2];
            if dat(1,1)+dat(1,2) > dat(3,1)+dat(3,2)
                ixY(3:4) = 1; 
            else
                ixY(3:4) = 3; 
            end
        elseif ixY(2)== 2 && ixX(2) == 3 % cases 2 & 4
            ixX(3:4) = [2,3];
            if dat(1,2)+dat(1,3) > dat(3,2)+dat(3,3)
                ixY(3:4) = 1; 
            else
                ixY(3:4) = 3; 
            end                            
        elseif ixX(2)== 2 && ixY(2) == 1 % cases 2 & 4
            ixY(3:4) = [1,2];
            if dat(1,1)+dat(2,1) > dat(1,3)+dat(2,3)
                ixX(3:4) = 1; 
            else
                ixX(3:4) = 3; 
            end
        elseif ixX(2)== 2 && ixY(2) == 3 % cases 2 & 4
            ixY(3:4) = [2,3];
            if dat(2,1)+dat(3,1) > dat(2,3)+dat(3,3)
                ixX(3:4) = 1; 
            else
                ixX(3:4) = 3; 
            end
        elseif abs(ixY(2)-2)+abs(ixX(2)-2) == 2  % cases 7 & 8
            ixX(3:4) = fliplr(ixX(1:2));
            ixY(3:4) = ixY(1:2);    
        end
        disp(iserr2by2_)
        iserr2by2_ = [];
    end
    linIx = sort(sub2ind(size(dat),ixY,ixX)); % linear index
    [ytemp1, xtemp1] = ind2sub([3 3],linIx(1));
    ytemp = pxYc + ytemp1 - 2;
    xtemp = pxXc + xtemp1 - 2;
    temp = zeros(m); temp(ytemp,xtemp)=1; temp = padarray(temp,[1 1]);
    Window2by2 = Window7by7(ytemp:ytemp+1, xtemp:xtemp+1);
    [yCoor2by2(i), xCoor2by2(i)] = find(temp==1); % rectangle coord for 2 by 2

    % FIT
%             numParam1_ = LB(:,i2) == 0; numParam2_ = UB(:,i2) == 0;
%             numParam_ = find(numParam1_ .* numParam2_ == 1);
    numParam_ = nSptsNonReg*5+1+5;

    if ~isempty(numParam_)
        numParam = numParam_(1);
        LB_ = LB(1:numParam,i2);
        UB_ = UB(1:numParam,i2);
        G0nonReg_ = G0nonReg(1:numParam,i2);
    else
        LB_ = LB(:,i2);
        UB_ = UB(:,i2);
        G0nonReg_ = G0nonReg(1:size(LB,1),i2);
    end
    tic
    [fitVal_,errFit7by7(i),err2res_,EXITFLAG] = gausFit(data7by7diff,funSpNonReg,G0nonReg_,LB_,UB_);
    time = toc;
    disp(sprintf('time:%.02f, #:%i, nSptsNonReg: %i,nSptsReg: %i, EXITFLAG: %i ',time,i,nSptsNonReg,nSptsReg,EXITFLAG));

    fitVal(:,i) = 0;
    fitVal(1:length(fitVal_),i) = fitVal_;

    %% calc XY value on the image (xyzDataGaus)
    Xgaus(i) =  double(xBW(i)) + fitVal(3,i) -px0;
    Ygaus(i) =  double(yBW(i)) + fitVal(5,i) -px0;
    INTnorm(i) = intPeak(i)*fitVal(2,i); % intensity

%            disp(sprintf('fr#:%i, i#:%i, window center:%.02f,%.02f, EXITFLAG: %i ',ixTrack(i),i,xBW(iyTrack(i),ixTrack(i)),yBW(iyTrack(i),ixTrack(i)),EXITFLAG(3)));

    %% for display only
    viewSpot2dispLoop;

    waitbar(i/nSpots);
    sp2(i) = i;
    i = i + 1;
end
