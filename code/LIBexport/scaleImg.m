% img imgMin,imgMax,imgGamma
nRow = nTileXY(1);
nCol = nTileXY(2);

[szY szX] = size(img);
sY = (szY-2*nRow)/nRow+2;
sX = (szX-2*nCol)/nCol+2;

slideValMin = get(hImgMin,'Value');    
slideValMax = get(hImgMax,'Value');    
minG = 0.2; maxG = 8.2; rG = maxG-minG;
slideValGamma = get(hImgGamma,'Value');    

for i = 1:nRow
    for j= 1:nCol
        img_ = img(sY*(i-1)+2:sY*(i-1)+sY-1,sX*(j-1)+2:sX*(j-1)+sX-1);
        imgMax_ = double(max(img_(:)));
        imgMin = slideValMin*imgMax_; 
        imgMax = slideValMax*imgMax_; 
        imgGamma  = slideValGamma*rG+minG;
        
        img_(find(img_<imgMin)) = imgMin;
        img_(find(img_>imgMax)) = imgMax;
        scale = double(max(img_(:)))^imgGamma/double(max(img_(:)));
        img_ = double(img_)/scale;
        img_ = img_.^imgGamma;
        img(sY*(i-1)+2:sY*(i-1)+sY-1,sX*(j-1)+2:sX*(j-1)+sX-1)=img_;
    end
end


