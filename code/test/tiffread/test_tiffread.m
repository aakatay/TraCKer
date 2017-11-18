fname = '2x2_200000frames.tif';
% fname = '128x88xy_1.9ms_32X25Y44x40_000.tif';
fname = '128x88xy_1.9ms_100P_300G_1.5x_000.tif';
t = Tiff(fname,'r');
t.getTag('ImageLength')
t.getTag('ImageWidth')

tiffDescrpt = t.getTag('ImageDescription')
ch1 = strfind(tiffDescrpt,'images');
ch2 = strfind(tiffDescrpt,'slices');
Frames = str2num(tiffDescrpt(ch1+7:ch2-1))



for i = 1:200000
    a(:,:,i)=tiffread
end


%% part2 create a tiff and read

