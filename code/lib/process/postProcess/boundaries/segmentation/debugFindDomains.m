
    load Rc1
    %[Rd, strInfo,Btb] = findDomains2(Rc1);
    [Rd, strInfo,Btb] = findDomains(Rc1);
    figure(1121);
    imagesc(Rd)
    mg=4;
    set(gca,'units','pixels','Position',[0 0 size(Rc,2)*mg size(Rc,1)*mg]); 
    set(gcf,'units','pixels','Position',[0 0 size(Rc,2)*mg size(Rc,1)*mg]); 
    imgFig = getframe(gcf);
    dataImg = imgFig.cdata; 
    imwrite(uint16(dataImg),'debugFindDomains.tif'); % structMap
    
    