BCKGRND = '120x132xy_10ms_100P_300G_1.5x_001_BACKGROUND.tif';


imageInfo=imfinfo(BCKGRND);
Frames=length(imageInfo);

for i = 1:Frames
    temp = imread(fname,i);
    bck(i) = median(median(temp));
end

plot(bck);
