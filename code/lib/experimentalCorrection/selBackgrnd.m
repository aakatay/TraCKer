%clear all;
%close;

%% select a roi for coefficient normalization
warning('off','MATLAB:class:EventWillBeRenamed')

load fname;
fTif = dir('*.tif');
for i = 1:length(fTif)
    if ~isempty(strfind(fTif(i).name,'MAX'))
        fname_ = [fTif(i).name];
    end
end
sprintf('background filename: %s',fname_)

%fname_ = ['MAX_' fname(1:end-8) '_000.tif'];
if isempty(dir(fname_))
    warning('no background file was found, no background correction will be implemented!!!')
    
    return;
    % a = input('do you want to skip intensity correction[y/n]:','s');
    if strcmp(a,'y')
        BWbckgrnd = ones;
    end
end
A = imread(fname_);
fnameMaxProj = fname_;
maxA = max(A(:));
A = maxA - A;
maxA = max(A(:));
Adisp = uint16(A-maxA/1.5);

figure('units','normalized','outerposition',[0 0 1 1]);
hImg = imagesc(Adisp);
colormap(gray)
axis image

hTextQuit = uicontrol('style','text','String','press the buttons and select regions','Position',[50 140 400 15]);
hBtnBck = uicontrol('style','pushbutton','units','pixel','position',[50 110 100 20],'Callback', 'q=1','String','select background'); 
hBtnROI = uicontrol('style','pushbutton','units','pixel','position',[180 110 60 20],'Callback', 'q=2','String','select ROI'); 
hBtnQUIT = uicontrol('style','pushbutton','units','pixel','position',[270 110 60 20],'Callback', 'q=3','String','QUIT'); 
hTextBckCoeff = uicontrol('style','text','String','quit','Position',[50 90 40 15]);
hBckCoeff = uicontrol('style','slider','units','pixel','position',[50 50 300 20]); 
%set(hBckCoeff,'Value',2); 

npol = 1;
q = 0;
while q ~= 3 % quit
    while q == 0
        listenSlider = addlistener(hBckCoeff,'ActionEvent',@(hObject, event) updBckCoeff(hObject, event,hImg,A,maxA,hTextBckCoeff)); 
        %listen2 = addlistener(hBtnQuit,'ActionEvent',@(hObject, event) updBckCoeff(hObject, event,hImg,A,maxA)); 
 
        pause(0.1)
    end
    if q == 1
        set(hTextQuit,'String','selecting background'); 
    elseif q == 2
        set(hTextQuit,'String','selecting ROI'); 
    elseif q == 3
        set(hTextQuit,'String','quitting'); 
    end
    if q ~= 3
        %while q == 0
        hPoly = imrect;
        pos = getPosition(hPoly);
        posRectCrop = pos;
        x = pos(1);
        y = pos(2);
        w = pos(3); 
        h = pos(4);
        posPoly = round([x y; x+w-1 y; x+w-1 y+h-1 ; x y+h-1]);
        delete(hPoly)
        hr=rectangle('Position', [x y w h],'FaceColor',[0 0.2 0.2 ],'LineStyle','none')
        q_=q; q=0;
    end
    if q_ == 1       
        BWbckgrnd(:,:,npol) = roipoly(ones(size(A)),posPoly(:,1),posPoly(:,2));
        posRectCropBCKGRND = posRectCrop;
    elseif q_ == 2
        BWroi = roipoly(ones(size(A)),posPoly(:,1),posPoly(:,2));;
        posRectCropROI = posRectCrop;
    end
end
    %npol = npol + 1;
%end
imagesc(sum(BWbckgrnd+BWroi,3)); axis image
save('BWselect','BWroi','BWbckgrnd','posRectCropBCKGRND','posRectCropROI','fnameMaxProj');

