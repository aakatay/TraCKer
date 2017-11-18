clear all; close all;
dbg=0;
load structMapUpd150421-cell10; ix = [6,8,12,18,20]; str='150421-cell10';
%load structMapUpd160726-cell8; ix = [7,15,19,20]; str='160726-cell8';

for i = 1: numel(ix)
    szXY = size(LtNew); szXY(3)=[];
    A = zeros(szXY);
    b0 = BtNew{ix(i)};
    [b]=extendBoundaries(b0,szXY);
    matFN  = [str '-struct' int2str(ix(i)) '.mat'];
    save(matFN,'b');
    if dbg
        
        imagesc(A);
        hold on
        plot(b0(:,2),b0(:,1),'g','LineWidth',1.5)
        plot(b(:,2),b(:,1),'r','LineWidth',1.5)
        hold off
        
    end
    cc=3;
end





        %A(sub2ind(szXY,b0(:,1),b0(:,2)))=1;
        %A(sub2ind(szXY,b(:,1),b(:,2)))=200;