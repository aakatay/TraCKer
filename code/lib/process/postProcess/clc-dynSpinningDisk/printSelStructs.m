
FNc = sprintf('structs\\%ic.tif',i);
FNd = sprintf('structs\\%id.tif',i);
if exist(FNc), delete(FNc); end;
if exist(FNd), delete(FNc); end;
for iii = 1:7
    imwrite(Asel(:,:,iii),FNc,'WriteMode','append','Compression', 'none') 
    imwrite(Dsel(:,:,iii),FNd,'WriteMode','append','Compression', 'none') 
end

FNprof = sprintf('structs\\%iintProf.tif',i);
figure(988)
imgFig = getframe(gcf); 
imgOut = imgFig.cdata;
imwrite(imgOut,FNprof,'Compression', 'none') 