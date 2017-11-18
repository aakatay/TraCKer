% prints all open figures

figHandles = findobj('Type','figure');
imgZFout3 = 'boundaryCalculation.tif';
delete(imgZFout3);

%sort
for i = 1: numel(figHandles)
    n(i)=figHandles(i).Number;
end
n=sort(n);
GCF=get(figure(n(2)));
sz= size(GCF.Children(1).Children.CData);
mag=4;
sz =sz*mag;
pos2= [0 0 sz(2)  sz(1)];
%print
for i = 1: numel(figHandles)-1

    figure(n(i));
colormap('gray')
    
    set(gcf,'units','pixels','Position',pos2); 
    set(gca,'units','pixels','Position',pos2); 

    imgFig = getframe(gcf); 
    imgOut = imgFig.cdata;
    imwrite(imgOut,imgZFout3,'WriteMode','append','Compression', 'none') 
end

    
colormap('parula')
    figure(n(end));
    
    set(gcf,'units','pixels','Position',pos2); 
    set(gca,'units','pixels','Position',pos2); 

    imgFig = getframe(gcf); 
    imgOut = imgFig.cdata;
    imwrite(imgOut,imgZFout3,'WriteMode','append','Compression', 'none') 