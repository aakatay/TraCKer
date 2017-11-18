function [Adev,Afamp,Amean] = findFlickering(fname,noiseThresh)
    % calculates the intensity dynamics of pre-bleach laptime movie
    infFN = imfinfo(fname);
    Frames = numel(infFN);
    for i = 1:Frames
        A(:,:,i) = imread(fname,i);
    end
    A = double(A);
    Adev =std(A,0,3); 
    Amean = mean(A,3);
    Afamp = max(A,[],3)-min(A,[],3); % flicker amplitude
    Adev(Afamp<noiseThresh) = 0; 
    
    % normalize
    mx = max([Amean(:)' Afamp(:)' Adev(:)']);
    Amean = Amean/mx;
    Afamp = Afamp/mx;
    Adev = Adev/mx;
    save mx mx
    
    fnameAdev = '..\lapDev.tif';
    fnameAmean = '..\lapMean.tif';
    fnameAfamp = '..\lapFamp.tif';
    imwrite(Adev,fnameAdev)
    imwrite(Amean,fnameAmean)
    imwrite(Afamp,fnameAfamp)
end
    
