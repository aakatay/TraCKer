function [ traceCoorDS ] = mapTraceCoor( traceCoor, rotAngle1, axis1, rotAngle2, axis2)
%MAPTRACEDATA calculates the x,y,z values based on the input rotation
% angles and axes
% traceCoor : trace coord. for 3Drotated images
% traceCoorDS : trace coord. for DS images

    X = traceCoor(:,1);
    Y = traceCoor(:,2);
    Z = traceCoor(:,3);
    
    
    
    traceCoorDS(1) = X;
    traceCoorDS(2) = Y;
    traceCoorDS(3) = Z;
    

end