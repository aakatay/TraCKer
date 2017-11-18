function imgC = bw2color(img,CG)
% CG :number of color grades
    
    CData = 1:CG;
    
    Zmin = min(CData(CData~=0));
    Zmax = max(CData(:));
    Zrange = Zmax - Zmin;
    Cplot = round((CG-1)*(CData-Zmin) / Zrange)+1;
    CM = colormap(jet(Zmax));    close;
    
    img = ceil(img/max(img(:))*(CG-1));
    imglin = img(:)+1;
    
    for i = 1:size(img,3)
         imgC_ = CM(imglin(:,i),:);
         imgC_ = reshape(imgC_,size(img,1),size(img,2),3);
         imgC(:,:,:,i) = imgC_;
    end
    plotColorBar(CData,5);  % ,'recruitment per 30 sec'
    fOut = 'blockMax6framesColor';
    colorbarPrint = sprintf('stats\\%s_colorbar.tiff',fOut);
    print(colorbarPrint,'-dtiff','-r80'); 

end