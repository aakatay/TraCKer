% extends to draw around pixels io from centers
% called inside findEdgeDistanceDistribution.m
function [b]=extendBoundaries(b,szXY)
    dbg = 0;
    
    yl = b(:,1);
    xl = b(:,2);
    
    %% increase border thickness in magnified resolution
    RbH = zeros(szXY*2);
    RbH0 = zeros(szXY);
    %xl = [floor(xl/2) xl];
    %yl = [floor(yl/2) yl];
    c=b*2;
    cw1=[1 0 1; 0 0 0; 1 0 1];
    cw2=[0 1 0 ; 1 0 1; 0 1 0];
    cw3=[1 0 0; 0 0 0; 0 0 1];
    cw4=[0 0 1; 0 0 0; 1 0 0];
    %plot(c(:,2),size(Rc,1)-c(:,1)+1,'LineWidth',1);
    c0x = 0;
    c0y = 0;
    for j= 1:size(c,1)
        RbH0(yl(j),xl(j))=1;
        RbH(c(j,1),c(j,2))=1;
    end

    if dbg, figure(916); imagesc(RbH0); axis image; end;

    if dbg, figure(917); imagesc(RbH); axis image; end;
    RbHc1=conv2(RbH,cw1,'same'); 
    RbHc2=conv2(RbH,cw2,'same');
    RbHc3=conv2(RbH,cw3,'same'); 
    RbHc4=conv2(RbH,cw4,'same');
    RbH(RbHc1==4)=1;
    RbH(RbHc2==2)=1;
    RbH(RbHc3==2)=1;
    RbH(RbHc4==2)=1;
    if dbg, figure(918); imagesc(RbH); axis image; end;
    RbH = conv2(RbH,ones(3),'same');
    RbH = im2bw(RbH,0);

    RbH = circshift(RbH,-[1 1]);
    
    if dbg, figure(919); imagesc(RbH); axis image; end;
    
    
    
    %% fill inside boundaries
    RbHws=watershed(RbH,4);
    if dbg, figure(920); imagesc(RbHws); axis image; end;

    y=[];
    RbHwsEl = unique(RbHws);
    for i = 1:length(RbHwsEl)
        if RbHwsEl(i)==1, continue; end
        y(i) = sum(sum(RbHws==RbHwsEl(i)));
    end
    [~,ix]=max(y);
    RbH(RbHws == RbHwsEl(ix))=1;
    
    %% find boundaries
    if dbg, figure(921); imagesc(RbH); axis image; end;
    [B,Lt]=bwboundaries(RbH,4);
    b=cell2mat(B(1));
    if dbg, figure(922); imagesc(RbH); axis image; 
        hold on;
        plot(b(:,2),b(:,1))
        hold off;
    end;
    RbH0mag=repelem(RbH0,2,2);
    
    %% display boundaries with corners
    b=cell2mat(B(1))+0.5;
    if dbg, figure(923); imagesc(RbH0mag); axis image; 
        hold on;
        plot(b(:,2),b(:,1))
        hold off;
    end;
    RbHdiff = RbH-RbH0mag;
    if dbg, figure(924); imagesc(RbHdiff); axis image; end;

    
    
    %% remove corners 
    cwin1 = ones(2);
    cwin2 = cwin1; cwin3 = cwin1; 
    cwin2(2) = 0; cwin3(3) = 0;
    cwin4 = -1*ones(4); cwin4(1:2,1:2)=1; cwin4(11)=1;
    cwin4 = fliplr(flipud(cwin4));
    cwin5 = -1*ones(3); cwin5(5)=1;
    
    % remove corners 1
    RbHcrnr = zeros(size(RbHdiff));
    RbHcrnrCv = conv2(RbHdiff,cwin1,'same')==4;
    if dbg, figure(925); imagesc(RbHcrnrCv); axis image; end;

    [B,Lt]=bwboundaries(RbHcrnrCv,4);

    y=[];
    RbHwsEl = unique(Lt);
    for i = 1:length(RbHwsEl)
        if RbHwsEl(i)==0, continue; end
        y(i) = sum(sum(Lt==RbHwsEl(i)));
    end
    [~,ix]=max(y);
    RbHcrnrCv(Lt == RbHwsEl(ix))=0;
    %if dbg, figure(999); imagesc(RbHcrnrCv); axis image; end;
    RbH = RbH-RbHcrnrCv;
    if dbg, figure(932); imagesc(RbH); axis image; end;
    % remove corners 2
    RbHcrnr = zeros(size(RbHdiff));
    RbHcrnrCv = conv2(RbHdiff,cwin2,'same')==3;
    if dbg, figure(933); imagesc(RbHcrnrCv); axis image; end;

    [B,Lt]=bwboundaries(RbHcrnrCv,4);

    y=[];
    RbHwsEl = unique(Lt);
    for i = 1:length(RbHwsEl)
        if RbHwsEl(i)==0, continue; end
        y(i) = sum(sum(Lt==RbHwsEl(i)));
    end
    [~,ix]=max(y);
    RbHcrnrCv(Lt == RbHwsEl(ix))=0;
    %if dbg, figure(999); imagesc(RbHcrnrCv); axis image; end;
    RbH = RbH-RbHcrnrCv;
    RbH(RbH<0)=0;
    if dbg, figure(934); imagesc(RbH); axis image; end;

    % remove corners 3
    RbHcrnr = zeros(size(RbHdiff));
    RbHcrnrCv = conv2(RbHdiff,cwin3,'same')==3;
    if dbg, figure(935); imagesc(RbHcrnrCv); axis image; end;

    [B,Lt]=bwboundaries(RbHcrnrCv,4);

    y=[];
    RbHwsEl = unique(Lt);
    for i = 1:length(RbHwsEl)
        if RbHwsEl(i)==0, continue; end
        y(i) = sum(sum(Lt==RbHwsEl(i)));
    end
    [~,ix]=max(y);
    RbHcrnrCv(Lt == RbHwsEl(ix))=0;
    %if dbg, figure(999); imagesc(RbHcrnrCv); axis image; end;
    RbH = RbH-RbHcrnrCv;
    RbH(RbH<0)=0;
    if dbg, figure(936); imagesc(RbH); axis image; end;

    % remove corners 4
    RbHcrnr = zeros(size(RbHdiff));
    if dbg, figure(937); imagesc(RbHdiff); axis image; end;

    RbHcrnrCv = conv2(RbHdiff,cwin4,'same')==5;
    RbHcrnrCv = circshift(RbHcrnrCv,[1 1]);
    if dbg, figure(938); imagesc(RbHcrnrCv); axis image; end;
    %if dbg, figure(925); imagesc(conv2(RbHdiff,cwin4,'same')); axis image; end;
    %if dbg, figure(999); imagesc(RbHcrnrCv); axis image; end;
    RbH = RbH-RbHcrnrCv;
    RbH(RbH<0)=0;
    if dbg, figure(932); imagesc(RbH); axis image; end;

    % remove corners 5
    RbHcrnr = zeros(size(RbHdiff));
    if dbg, figure(939); imagesc(RbHdiff); axis image; end;

    RbHcrnrCv = conv2(RbHdiff,cwin5,'same')==1;
    %RbHcrnrCv = circshift(RbHcrnrCv,[1 1]);
    if dbg, figure(925); imagesc(RbHcrnrCv); axis image; end;
    %if dbg, figure(925); imagesc(conv2(RbHdiff,cwin4,'same')); axis image; end;
    %if dbg, figure(926); imagesc(RbHcrnrCv); axis image; end;
    RbH = RbH-RbHcrnrCv;
    RbH(RbH<0)=0;
    if dbg, figure(940); imagesc(RbH); axis image; end;
        
    
    
    %% display updated 
    [B,Lt]=bwboundaries(RbH,4);
    if dbg, figure(1001); imagesc(RbH); axis image; 
        hold on;
        b=cell2mat(B(1));
        plot(b(:,2),b(:,1))
        hold off;
    end;
    RbH0mag=repelem(RbH0,2,2);
    
    b=cell2mat(B(1))+0.5;
    if dbg, figure(1002); imagesc(RbH0mag); axis image; 
        hold on;
        plot(b(:,2),b(:,1),'r')
        hold off;
    end;
    b=b/2+0.25;
end