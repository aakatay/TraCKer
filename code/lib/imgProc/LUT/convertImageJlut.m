ijout = 'ijout.tif';
if ~exist(ijout)
    a= 1:256;
    imwrite(uint16(a),'ijIN.tif');
else
    b= imread(ijout);
end


