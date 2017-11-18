% calculates the positions based on the set of the defined points
function [posIx, pos] = cellGetOutline_calcPos(posIx,pos)
    Npos = numel(posIx);
    Np = size(pos,1); % # points
    if posIx(end) ~= size(pos,3)
        pos(:,:,end) = pos(:,:,posIx(end));
        posIx(end) = size(pos,3);
    end
    for i = 1:Npos-1    
        if posIx(i+1)-posIx(i) > 1
            dN = (posIx(i+1)-posIx(i));
            dX = (pos(:,1, posIx(i+1))-pos(:,1,posIx(i))) / dN;
            dY = (pos(:,2,posIx(i+1))-pos(:,2,posIx(i))) / dN;
            m = 1;
            for j = posIx(i)+1: posIx(i+1)
                pos(:,1,j) = pos(:,1,posIx(i))+dX*m;
                pos(:,2,j) = pos(:,2,posIx(i))+dY*m; 
                m=m+1;
            end
        end
    end
end