function setTiffDescription(varargin)
    %return
    tifname = cell2mat(varargin(1));
    tifDescript = cell2mat(varargin(2));
    isModifiedImage=0;
    if nargin == 3
        isModifiedImage = cell2mat(varargin(3));
    end
    
    if isModifiedImage
        tifDescript = ['IMAGE IS PROCESSED (see below). Original Image Info:' tifDescript];
    end
    t = Tiff(tifname,'r+');
    t.setTag('ImageDescription',tifDescript);
    rewriteDirectory(t)
    t.close();
end


