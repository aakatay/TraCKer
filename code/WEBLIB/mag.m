function [magnitude] = mag(vect)
% calculates the magnitude of the vector
    vect = abs(vect);
    if size(vect,1)+size(vect,2)-1 == numel(vect)
        if numel(vect) == 3
            magnitude = sqrt(vect(1)^2 + vect(2)*vect(2) + vect(3)^2);
        elseif numel(vect) == 2
            magnitude = sqrt(vect(1)^2 + vect(2)^2);
        end
    else
        if size(vect,2) == 3
            magnitude = sqrt(vect(:,1).^2 + vect(:,2).^2 + vect(:,3).^2);
        elseif size(vect,2) == 2
            magnitude = sqrt(vect(:,1).^2 + vect(:,2).^2);
        else
            error('mag.m: incorrect input size')
        end
    end
end