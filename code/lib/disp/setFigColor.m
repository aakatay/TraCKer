function setFigColor(hFig,hAx)
%return
set(hFig,'Color',[0 0 0]); 
typ = get(hAx(1),'Type');
if strcmp(typ,'line') || strcmp(typ,'axes')
    set(hAx(1),'Color',[0 0 0]);
end
set(hAx(1),'YColor',[1 1 1],'XColor',[1 1 1]);

if numel(hAx)==2
    set(hAx(2),'XColor',[1 1 1]);
end