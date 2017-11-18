ix = [12,18,24,41,54,58,65,115,125,131]; % selected indices
ix = [366,691,2846,2849]; % selected indices
ix = 15754;

FNADconc = 'ADconc.tif';
FNADdiff = 'ADdiff.tif';
if exist(FNADconc), delete(FNADconc); end
if exist(FNADdiff), delete(FNADdiff); end
for j = 1:numel(ix)
    i=ix(j);
    FNc = sprintf('structs\\%ic.tif',i);
    FNd = sprintf('structs\\%id.tif',i);
    FNprof = sprintf('structs\\%iintProf.tif',i);
    
    FNcConc = sprintf('structs\\%icConc.tif',i);
    FNdConc = sprintf('structs\\%idConc.tif',i);
    
    
    sp = ones(9,1)*256;
    Aconc=sp;
    Dconc=sp;
    for iii = 1:7
        Asel(:,:,iii)= imread(FNc,iii);
        Dsel(:,:,iii)= imread(FNd,iii); 
        Aconc = [Aconc Asel(:,:,iii) sp];
        Dconc = [Dconc Dsel(:,:,iii) sp];
    end
    sp2 = ones(3,size(Aconc,2))*256;
    ADconc = [sp2; Aconc ; sp2; Dconc;sp2];
    
    Amean1 = mean(Asel(:,:,1:3),3);
    Amean2 = mean(Asel(:,:,5:7),3);
    Adiff = Amean1 - Amean2;
    Dmean1 = mean(Dsel(:,:,1:3),3); 
    Dmean2 = mean(Dsel(:,:,5:7),3);
    Ddiff = Dmean1 - Dmean2;
    
    
    Anorm = double(Amean1);
    mx = max(Anorm(:));
    mn = min(Anorm(:));
    Anorm = (Anorm-mn)/(mx-mn)*256;
    Amean1=Anorm;
    
    Anorm = double(Amean2);
    mx = max(Anorm(:));
    mn = min(Anorm(:));
    Anorm = (Anorm-mn)/(mx-mn)*256;
    Amean2=Anorm;
    
    Anorm = double(Adiff);
    mx = max(Anorm(:));
    mn = min(Anorm(:));
    Anorm = (Anorm-mn)/(mx-mn)*256;
    Adiff=Anorm;
    
    sp3 = ones(size(Amean1,1),4)*256;
    Ad = [sp3 Amean1 sp3 Amean2 sp3 Adiff sp3];
    Dd = [sp3 Dmean1 sp3 Dmean2 sp3 Ddiff sp3];
    ADd = [ones(1,size(Ad,2))*256; Ad; ones(3,size(Ad,2))*256; Dd; ones(1,size(Ad,2))*256];
    
    imwrite(ADconc,FNADconc,'WriteMode','append','Compression', 'none') 
    imwrite(uint16(ADd),FNADdiff,'WriteMode','append','Compression', 'none') 
    
    imwrite(Aconc,FNcConc);
    imwrite(Dconc,FNdConc);
    
end