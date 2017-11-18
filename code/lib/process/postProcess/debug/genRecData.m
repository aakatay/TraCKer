
A0m = repelem(A0,4,4);
A0m=circshift(A0m,sxsy(1),2);
Ad=circshift(A0m,sxsy(2),1);
fnameSum = 'binImgRcrtSum_genRecData.tif';
imwrite(Ad,fnameSum);
