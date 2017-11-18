% run twice
clear
    a = zeros(20,20);
    sh = 2;
    sh2 = 0;
    a(7+sh,10,1)=1; a(5+sh,10,1)=1; a(6+sh,8,1)=1;
    a(4+sh2,13,1)=1;a(6+sh2,13,1)=1;
    a(13,6,1)=1;a(12,8,1)=1;

    a(5+sh2,6,1)=1;a(7+sh2,6,1)=1;
    a(13,11,1)=1;a(12,13,1)=1;a(10,13,1)=1;

     
     
     if exist('inputSpots.mat')
         load('inputSpots');
         x = PxBW_;
         y = PyBW_;
         a = zeros(128);
         numSpots = size(x,2);
         frstSpot = 1;
         numSpots = 15;
         frstSpot = 10;
         numSpots = 11;
         spots = frstSpot:numSpots;
         for i = spots
             a(y(i),x(i)) = 1;
         end
     end
     
    a = a*1000;
    BWbckgrnd = zeros(size(a));
    BWbckgrnd(7:13,7:13)=1;
    save('BWbckgrnd','BWbckgrnd')
    bckImg = ones(size(a));
    imwrite(uint16(a),'regressSpotDetect_10ms_001.tif')
    imwrite(uint16(bckImg),'regressSpotDetect_10ms_000.tif')
    save_fname(1);
    try 
        TraCKer_3D_w_ZcolorPlot_deep
    catch err
        
        overlayIMGandDET
        maximize
        return
    end


    try 
        TraCKer_3D_w_ZcolorPlot_deep

    catch    
        ;
    end

    overlayIMGandDET
    maximize
    

    

