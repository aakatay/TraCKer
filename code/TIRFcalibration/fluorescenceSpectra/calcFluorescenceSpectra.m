% called by dispFluorescenceSpectra.m
%% calculates transmission of the spectra for each channel
% m : emission
% x : excitation
% FX : filter excitation
% FD : filter dichroic
% FDq : filter dichroic quadview
% FM : filter emission
% xFMq: filter emission quadview

dbg = 0;

% output
fnChanInt = 'channelIntensity';
if isBRT
    fnChanInt = ['BRT-' fnChanInt];
else
    fnChanInt = ['PRB-' fnChanInt];
end
   
if isCy3
    fnChanInt = ['CY3-' fnChanInt];
else
    fnChanInt = ['GFP-' fnChanInt];
end 

n = numel(xq);
P = ones(n,1); % initial intensity
%% LASER SOURCE
lw = 4; %laser band
for i = 1:laserN % each laser
    L = xq*0;
    [~,ixmin] = min(abs(xq-laserX(i)));
    L(ixmin) = 1;
    L = conv(L,ones(lw,1),'same');
    LS(:,i) = L/sum(L(:)); % normalize
end

if dbg
    dbgOffset = repmat(1:laserN,n,1);
    figure(31); plot(xq,dbgOffset+LS*0.9);
end

%% FLUO
% excitation
for i = 1:fluoN % each fluorophore
    FLx(:,i) = flDx{i};
end

% emission
for i = 1:fluoN % each fluorophore
    FLm(:,i) = flDm{i};
end

%% (chls & chfl)  EXCITATION & EMISSION TRANSFER
if isBRT % quadview
    % (L) Laser transmission
    CH1X = [P      fcD{1}      1-fcD{3}    fcD{2}    1-fqd.D{1}        fqd.D{2}*0+1    fqd.D{3}*0+1    fqm.D{1}        fqm.D{2}*0+1    fqm.D{3}*0+1]; % CH1
    CH2X = [P      fcD{1}      1-fcD{3}    fcD{2}      fqd.D{1}      1-fqd.D{2}        fqd.D{3}*0+1    fqm.D{1}*0+1    fqm.D{2}        fqm.D{3}*0+1]; % CH2
    CH3X = [P      fcD{1}      1-fcD{3}    fcD{2}      fqd.D{1}        fqd.D{2}      1-fqd.D{3}        fqm.D{1}*0+1    fqm.D{2}*0+1    fqm.D{3}]; % CH3
          %[x0     xFX           xFD       XFM         xFDq1           xFDq2           xFDq3           xFMq1           xFMq2           xFMq3];
    % (F) fluo emission transmission
    CH1M = [P      fcD{1}*0+1    fcD{3}    fcD{2}    1-fqd.D{1}        fqd.D{2}*0+1    fqd.D{3}*0+1    fqm.D{1}        fqm.D{2}*0+1    fqm.D{3}*0+1]; % CH1
    CH2M = [P      fcD{1}*0+1    fcD{3}    fcD{2}      fqd.D{1}      1-fqd.D{2}        fqd.D{3}*0+1    fqm.D{1}*0+1    fqm.D{2}        fqm.D{3}*0+1]; % CH2
    CH3M = [P      fcD{1}*0+1    fcD{3}    fcD{2}      fqd.D{1}        fqd.D{2}      1-fqd.D{3}        fqm.D{1}*0+1    fqm.D{2}*0+1    fqm.D{3}]; % CH3
          %[m0     xFX           xFD       XFM         xFDq1           xFDq2           xFDq3           xFMq1           xFMq2           xFMq3];

    ch1x = cumprod(CH1X,2);
    ch2x = cumprod(CH2X,2);
    ch3x = cumprod(CH3X,2);
    ch1m = cumprod(CH1M,2);
    ch2m = cumprod(CH2M,2);
    ch3m = cumprod(CH3M,2);

    if dbg
        np = size(CH2X,2); % num param
        dbgOffset = repmat(1:np,n,1);
        figure(11); plot(xq,dbgOffset+CH2X*0.9);
        hold on;
        plot(xq,dbgOffset+ch1x*0.9,'LineWidth',1,'Color',[0 0 0]);
        hold off;    
        ylim([1,np+1])
        figure(12); plot(xq,dbgOffset+CH2M*0.9);
        hold on;
        plot(xq,dbgOffset+ch2m*0.9,'LineWidth',1,'Color',[0 0 0]);
        hold off;    
        ylim([1,np+1])
    end

    CHL(:,1,1) = ch1x(:,end).*LS(:,1); % channel1 laser1
    CHL(:,1,2) = ch1x(:,end).*LS(:,2);
    CHL(:,1,3) = ch1x(:,end).*LS(:,3);
    CHL(:,2,1) = ch2x(:,end).*LS(:,1);
    CHL(:,2,2) = ch2x(:,end).*LS(:,2);
    CHL(:,2,3) = ch2x(:,end).*LS(:,3);
    CHL(:,3,1) = ch3x(:,end).*LS(:,1);
    CHL(:,3,2) = ch3x(:,end).*LS(:,2);
    CHL(:,3,3) = ch3x(:,end).*LS(:,3);

    CHF(:,1,1) = ch1m(:,end).*FLm(:,1); % channel1 fluo1
    CHF(:,1,2) = ch1m(:,end).*FLm(:,2);
    CHF(:,1,3) = ch1m(:,end).*FLm(:,3);
    CHF(:,2,1) = ch2m(:,end).*FLm(:,1);
    CHF(:,2,2) = ch2m(:,end).*FLm(:,2);
    CHF(:,2,3) = ch2m(:,end).*FLm(:,3);
    CHF(:,3,1) = ch3m(:,end).*FLm(:,1);
    CHF(:,3,2) = ch3m(:,end).*FLm(:,2);
    CHF(:,3,3) = ch3m(:,end).*FLm(:,3);

    CHL2D = reshape(CHL,[n, 9]);
    chls = reshape(sum(CHL2D,1),[3,3]); %[ch x ls]
    CHF2D = reshape(CHF,[n, 9]);
    chfl = reshape(sum(CHF2D,1),[3,3]); %[ch x fl]
else % PRB: CSU    (only one view channel)
    % (L) Laser transmission
    CH1X = [P      fcD{1}      1-fcD{3}    fcD{2}    fqm.D{1}]; % CH1
          %[x0     xFX           xFD       XFM       xFMq1   ];
    % (F) fluo emission transmission
    CH1M = [P      fcD{1}*0+1    fcD{3}    fcD{2}    fqm.D{1}]; % CH1
          %[m0     xFX           xFD       XFM       xFMq1   ];    
          
    ch1x = cumprod(CH1X,2);
    ch1m = cumprod(CH1M,2);

    CHL(:,1,1) = ch1x(:,end).*LS(:,1); % channel1 laser1
    CHL(:,1,2) = ch1x(:,end).*LS(:,2);
    CHL(:,1,3) = ch1x(:,end).*LS(:,3);

    CHF(:,1,1) = ch1m(:,end).*FLm(:,1); % channel1 fluo1
    CHF(:,1,2) = ch1m(:,end).*FLm(:,2);
    CHF(:,1,3) = ch1m(:,end).*FLm(:,3);
    
    chls = sum(CHL,1); %[ch x ls]
    chfl = sum(CHF,1); %[ch x fl]
    
end


if dbg
    nCH = size(CHL2D,2);
    dbgOffset = repmat(1:nCH,n,1);
    figure(21); plot(xq,dbgOffset+CHL2D*0.9);
    hold on;
    plot(xq,dbgOffset+CHF2D*0.9);
    hold off;
    ylim([1,nCH+1])
end

%% (flx) FLmUO EXCITATION 
for i = 1:laserN % each laser
    for j = 1:fluoN
        flx(i,j) = sum(LS(:,i).*FLx(:,j).*fcD{1});
    end
end


save(fnChanInt,'chls','chfl','flx')




