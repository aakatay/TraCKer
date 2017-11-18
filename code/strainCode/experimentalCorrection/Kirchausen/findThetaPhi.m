function  phiTheta = findThetaPhi(imgBase,xyzS)
% based on xyz coor of the input finds the closest spot and calc tilt angle

    filtMtrx = 0*imgBase;
    if exist('baseSpots.mat') && xyzS(1) == 0 
        load('baseSpots.mat');
        display('baseSpots.mat is loaded');
    else
        s = [5 5 5]; % spot size in px 
        sC = (s-1)/2; % for cropping after conv
        convMtrx = ones(s(1),s(2),s(3));
        for iS = 1:size(xyzS,2)
            xyzSMtrx = 0*imgBase;
            xyzSMtrx(xyzS(2,iS),xyzS(1,iS),xyzS(3,iS)) = 1;
            filtMtrx = convn(xyzSMtrx,convMtrx);
            filtMtrx = filtMtrx(sC(1)+1:end-sC(1),sC(2)+1:end-sC(2),sC(3)+1:end-sC(3)); % resize by cropping after convolution
            imgBase2 = filtMtrx.*imgBase;
            %figure(11); imagesc(imgBase2(:,:,xyzS(3,iS)));            
            [maxval, maxind] = max(imgBase2(:));
            [X(iS), Y(iS), Z(iS)] = ind2sub(size(imgBase2),maxind);


        end
        save('baseSpots.mat','X','Y','Z');
        devS=mean(abs(xyzS(3,:)-Z)+abs(xyzS(2,:)-X)+abs(xyzS(1,:)-Y)); % deviation of the selected center wrt to center of the spot
        if devS > 1
            sprintf('select better spots: %d',devS);
        end
    end
    [n V p] = affine_fit([X' Y' Z']);
    display('normal'); n
    
    phi     = atan(n(1)/n(2))/pi*180;
   
    theta   = acos(n(3))/pi*180*sign(phi);
    phiTheta = [phi theta];
end


    %mxPj = maxPventral/max(maxPventral(:));
    %nSpot   = 20;
    %Base = findLateralSpots(mxPj,MxZeInd,nSpot); % base of the cell 
    