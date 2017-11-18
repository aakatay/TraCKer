%Written by:    Joshua Ferguson
%               Kural Group
%               OSU Physics
%Last update:   8/6/13
function [IMG,DATA] = Gaussian_TIFF(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main Function: Will send relavant variable to workspace.
%filname is the filename of the tiff file, including extension
%frames is a 1D array containing the starting a stopping frames
%IMG is a 3D array. The first two dimensions is the intensity at given
%(x,y) coordinates. The third dimension is the frame from the TIFF movie
%that it is taken from.
%DATA is a 2D array where the second dimention coresponds to frame and 
%DATA(1,:) is x coordinate
%DATA(2,:) is y coordinate
%DATA(3,:) is x standard deviation
%DATA(4,:) is y standard deviation
%DATA(5,:) is amplitude.
%DATA(6,:) is the background, defined at the top of Gaussian2DJILLA
%DATA(7,:) is highest intensity pixel minus the background.
% Usage example:
% >> [IMG,DATA] = Gaussian_TIFF('Clathrin_PSF.tif',75);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 1
        IMG = cell2mat(varargin(1));
        frames= size(IMG,3);
    elseif nargin == 2
        filename = cell2mat(varargin(1));
        frames = cell2mat(varargin(2));
        for i=1:frames
        IMG(:,:,i) = imread(filename,i);
        end
    end

    
    for i=1:frames
    IMG = double(IMG);
    [DATA(1,i),DATA(2,i),DATA(3,i),DATA(4,i),DATA(5,i),DATA(6,i),DATA(7,i)] = Gaussian2DJILA(IMG(:,:,i));
    end
end

function [cx,cy,sx,sy,PeakOD,bkgd,pOD] = Gaussian2DJILA(m)
% following is adapted from JILA - but at this point quite far removed
% i.e. uses a different fitting function and more updated matlab syntax
% http://jilawww.colorado.edu/bec/BEC_for_everyone/matlabfitting.htm
% a function to fit a 2-D monochrome fluorescent TIFF image
% possible update: allow for colored images.
        
        options = optimset('Display','off');
    
        % find the center of mass coordinates and peak
        [sizey, sizex] = size(m);
        bkgd = min([mean(m,1),mean(m,2)']);
        m = m - bkgd;
        [cx,cy,sx,sy] = centerofmass(m);
        pOD = max(max(m));
        mOD = min(min(m));

        % estimate the 1D "x" coordinate
        mx = double(m(round(cy),:));
        x1D = double(1:sizex);
        ip1D = [double(cx),double(sx),double(pOD)];
        lb = [1,0,double(mOD)];
        ub = [double(sizex),inf,double(pOD)];
        fp1D = lsqcurvefit(@gaussian1D,ip1D,x1D,mx,lb,ub,options);
        cx = fp1D(1);
        sx = fp1D(2);
        PeakOD = fp1D(3);
        
        % estimate the 1D "y" coordinate
        my = double(m(:,round(cx))');
        y1D = double(1:sizey);
        ip1D = [double(cy),double(sy),double(PeakOD)];
        lb = [1,0,double(mOD)];
        ub = [double(sizey),inf,double(pOD)];
        fp1D = lsqcurvefit(@gaussian1D,ip1D,y1D,my,lb,ub,options);
        cy = fp1D(1);
        sy = fp1D(2);
        PeakOD = fp1D(3);
     
        % now perform 2D fit based on estimates of 1D fittings
        rcx = round(cx);
        rcy = round(cy);
        hwin= 2.;
        win = 2.*hwin + 1.;
        
        initpar = [double(cx),double(cy),double(sx),double(sy),double(PeakOD)];
        % getting a good window for evaluation
        if (hwin<rcx && rcx<(sizex-hwin) && hwin<rcy && rcy<(sizey-hwin))
            m = wkeep(m,[win win],[(rcy-hwin) (rcx-hwin)]);
            [XY(:,:,1),XY(:,:,2)] = meshgrid(double(rcx-hwin:rcx+hwin),double(rcy-hwin:rcy+hwin));
        elseif (hwin<rcy && rcy<(sizey-hwin) && rcx>(sizex-(hwin+1.)))
            m = wkeep(m,[win win],[(rcy-hwin) (sizex-(win-1.))]);
            [XY(:,:,1),XY(:,:,2)] = meshgrid(double(sizex-(win-1.):sizex),double(rcy-hwin:rcy+hwin));
        elseif (hwin<rcx && rcx<(sizex-hwin) && rcy>(sizey-(hwin+1.)))
            m = wkeep(m,[win win],[(sizey-(win-1.)) (rcx-hwin)]);
            [XY(:,:,1),XY(:,:,2)] = meshgrid(double(rcx-hwin:rcx+hwin),double(sizey-(win-1.):sizey));
        elseif (hwin<rcy && rcy<(sizey-hwin) && rcx<(hwin+1.))
            m = wkeep(m,[win win],[(rcy-hwin) 1.]);
            [XY(:,:,1),XY(:,:,2)] = meshgrid(double(1:win),double(rcy-hwin:rcy+hwin));
        elseif (hwin<rcx && rcx<(sizex-hwin) && rcy<(hwin+1.))
            m = wkeep(m,[win win],[1. (rcx-hwin)]);
            [XY(:,:,1),XY(:,:,2)] = meshgrid(double(rcx-hwin:rcx+hwin),double(1:win));
        elseif (rcx<(hwin+1.) && rcy<(hwin+1.))
            m = wkeep(m,[win win],[1. 1.]);
            [XY(:,:,1),XY(:,:,2)] = meshgrid(double(1:win),double(1:win));
        elseif (rcx>(sizex-(hwin+1.)) && rcy>(sizey-(hwin+1.)))
            m = wkeep(m,[win win],[(sizey-(win-1.)) (sizex-(win-1.))]);
            [XY(:,:,1),XY(:,:,2)] = meshgrid(double(sizex-(win-1.):sizex),double(sizey-(win-1.):sizey));
        elseif (rcy>(sizey-(hwin+1.)) && rcx<(hwin+1.))
            m = wkeep(m,[win win],[(sizey-(win-1.)) 1.]);
            [XY(:,:,1),XY(:,:,2)] = meshgrid(double(1:win),double(sizey-(win-1.):sizey));
        elseif (rcy<(hwin+1.) && rcx>(sizex-(hwin+1.)))
            m = wkeep(m,[win win],[1. (sizex-(win-1.))]);
            [XY(:,:,1),XY(:,:,2)] = meshgrid(double(sizex-(win-1.):sizex),double(1:win));
        end
        m = double(m);  % important background subtraction
        lb = [1,1,0,0,0];
        ub = [double(sizex),double(sizey),inf,inf,double(pOD)];
        fp = lsqcurvefit(@gaussian2D,initpar,XY,m,lb,ub,options);
        cx = fp(1);
        cy = fp(2);
        sx = fp(3);
        sy = fp(4);
        PeakOD = fp(5);
        

end

function [cx,cy,sx,sy] = centerofmass(m)
% PURPOSE: find c of m of distribution

    [sizey, sizex] = size(m);
    vx = sum(m,1);
    vy = sum(m,2)';
    
    x = 1:sizex;
    y = 1:sizey;

    cx = sum(vx.*x)/sum(vx);
    cy = sum(vy.*y)/sum(vy);

    sx = sqrt(sum(vx.*((x-cx).^2))/sum(vx));
    sy = sqrt(sum(vy.*((y-cy).^2))/sum(vy));
end

function [z] = gaussian1D (p,x)
%p(1) is center of mass
%p(2) is standard deviation
%p(3) is max value
z = p(3)*exp(-0.5*(((p(1)-x)/p(2)).^2));
%z = sum(zx.^2);
end

function [z] = gaussian2D(p,XY)
%all below have been adapted by previous fittings
%p(1) is center of mass of x;
%p(2) is standard deviation of x;
%p(3) is center of mass of y;
%p(4) is standard deviation of y;
%p(5) is the peak value;

    z = p(5)*exp(-0.5*(((XY(:,:,1)-p(1))/p(3)).^2+((XY(:,:,2)-p(2))/p(4)).^2));
end
