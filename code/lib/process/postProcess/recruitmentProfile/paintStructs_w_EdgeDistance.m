
stFN = rdir('structMapUpd*.mat');
edgeFN = rdir('edgeDistData*.mat');
psFN = 'edgeDistPaintStructs.tif';
psFigFN = 'edgeDistPaintStructs.fig';

if isempty(edgeFN), return; end;
load(stFN.name)
load(edgeFN.name)

r = edgeDistMeanArr(:,1)./edgeDistMeanArr(:,2);

rbck = min(r)-0.1;
ns = size(r,1);
for i=1:ns
   LtNew(:,:,i) = LtNew(:,:,i)*r(i);
end
Lt = sum(LtNew,3);
Lt(Lt==0)=rbck;
sx = size(Lt,2);
sy = size(Lt,1);
imagesc(Lt);
colorbar; axis image
set(gcf,'units','pixels','Position',[200,120,sx+100,sy]); 
set(gca,'units','pixels','Position',[0,0,sx,sy]);
CM=colormap('parula');
CM(1,:)=0;
colormap(CM)

%delete(psFN)
if ~exist('paintStructs_w_EdgeDistance_dontPrint')
    imgFig = getframe(gcf);
    dataImg = imgFig.cdata; 
    imwrite(uint16(dataImg),psFN,'Compression', 'none'); % structMap
end