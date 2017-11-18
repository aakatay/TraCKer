clear
close all

fROI = fopen('recComp\ROI.txt');

fcoor = [];
%% check the output files
fdata = rdir('recComp\\data*.mat');
fcoor_ix = [];
for i = 1:numel(fdata)
    fn = fdata(i).name;
    ixDot = strfind(fn,'_');
    ixDot = ixDot(1);
    fcoor{i} = fn(13:ixDot-1); % file index
    fcoor_ix(i) = str2num(fn(end-5:end-4));
end
fnameSum = rdir('binImgRcrtSum*.tif');
fnameSum = fnameSum(2).name; % run recruitComp first

A=double(imread(fnameSum));
A = double(im2bw(A));
figure(1);imagesc(A); axis image; 
figure(2);imagesc(zeros(size(A))); axis image; 

s=0;n=0;p=0;t=0;h=0;j=0;m=0;ns=0;
while 1
    Col = 'k'; % color
    coors = fgetl(fROI);
    coorsAr{s+1} = coors;
    if sum(coors == -1) || isempty(coors), break, end % eof
     isNonStruct =0;
    if strncmp(coors,'p',1)
        coors = coors(2:end);
        isPit = 1; isNonStruct =1;
        p=p+1;
        Col = 'g'; 
    elseif strncmp(coors,'n',1) % new structures
        coors = coors(2:end);
        isNew = 1; isNonStruct =1;
        n=n+1; 
        Col = 'y'; 
    elseif strncmp(coors,'c',1) % joining structures
        coors = coors(2:end);
        isJoin = 1; isNonStruct =1;
        j=j+1; 
        Col = 'c'; 
    elseif strncmp(coors,'t',1) % tearing structures
        coors = coors(2:end);
        isTear = 1; isNonStruct =1;
        t=t+1; 
        Col = 'r'; 
    elseif strncmp(coors,'h',1) % hot spots
        coors = coors(2:end);
        isHot = 1; isNonStruct =1;
        h=h+1; 
        Col = 'w'; 
    elseif strncmp(coors,'m',1) % streaming
        coors = coors(2:end);
        isStream = 1; isNonStruct =1;
        m=m+1; 
        Col = 'r'; 
    elseif strncmp(coors,'g',1) % gripping
        coors = coors(2:end);
        isStream = 1; isNonStruct =1;
        m=m+1; 
        Col = 'm'; 
    end
    a=sscanf(coors,'%dX%dY%dx%d');
    disp(coors)
    fix = find(strcmp(fcoor,coors));
    fix_ix = fcoor_ix(fix);
    if isNonStruct ns = ns + 1; end
    if isempty(fix) &&  ~isNonStruct;
        warning(sprintf('ROI not processed: data%s.mat',coors))
    elseif isNonStruct
        ;
    else
        ixStay = find(~strcmp(fcoor,coors));
        fcoor = fcoor(ixStay);
        fcoor_ix = fcoor_ix(ixStay);
    end
    xCr = a(1)*4; yCr = a(2)*4;  szX = a(3)*4; szY = a(4)*4;
    A(yCr+1:yCr+szY,xCr+1:xCr+szX) = 0; % structure image
    
    figure(2) % showROIix indices
    rectangle('Position',[xCr yCr szX szY],'EdgeColor',Col)
    tx=xCr+szX/2;
    ty=yCr+szY/2;
    th=text(tx-5,ty,sprintf('%02i',fix_ix));
    set(th,'Color',[1 0 0]);
    
    figure(1) %showROIix 
    rectangle('Position',[xCr yCr szX szY],'EdgeColor',Col)
    
    s=s+1;
    
end
if ~isempty(fcoor)
    warning('delete the following ROI data:')
    fcoor
end

disp(sprintf('# of structures:%i, # of non-structures:%i',s-ns,ns))
imwrite(uint16(A),'binImgRcrtSum_ROIdeleted.tif')
%% print showROI
figure(1)
[szY,szX]=size(A);
set(gcf,'units','pixels','Position',[200,200,szX+50,szY+120]); 
set(gca,'units','pixels','Position',[50,50,szX,szY]); axis tight;
title('structures:: pits(p):green, new(n):yellow, joining(c):cyan, tearing(t):red, hotspots(h):white, streaming(m):red, gripping(g):magenta stable:black')
imgFig = getframe(gcf);
imgData = imgFig.cdata; 
imwrite(uint16(imgData),'recComp\showROI.tif')

%% print showROIix indices
figure(2)
set(gcf,'units','pixels','Position',[200,200,szX+50,szY+120]); 
set(gca,'units','pixels','Position',[50,50,szX,szY]); axis tight;
title('structures:: pits(p):green, new(n):yellow, joining(c):cyan, tearing(t):red, hotspots(h):white, streaming(m):red, gripping(g):magenta stable:black')
imgFig = getframe(gcf);
imgData = imgFig.cdata; 
imwrite(uint16(imgData),'recComp\showROIix.tif')


