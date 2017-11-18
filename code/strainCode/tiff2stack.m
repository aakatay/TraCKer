% reads sequences of Tiff file and write RGB tiff stack
% genName : input generic name 
% fOut : output file

function tiff2stack(varargin) % genName,fOut
    
    if nargin <1
        genName = 'detectedSpots';
        genName = [];
        fOut = 'fileOut_tiff2stack.tif';
    elseif nargin == 1
        genName = cell2mat(varargin(1));
        fOut = [genName '_stack.tif'];
    else
        genName = cell2mat(varargin(1));
        fOut = cell2mat(varargin(2));
        
    end
        
    fname = dir([genName '*.tif']);
    nF = numel(fname); 
    
    for ixFrm = 1:nF        
        img3D = imread(fname(ixFrm).name,1);
        if ixFrm == 1
            imwrite(uint16(img3D),char(fOut),'tiff','Compression', 'none' );
        else
            imwrite(uint16(img3D),char(fOut),'tiff','WriteMode','append','Compression', 'none' );
        end
        delete(fname(ixFrm).name)
    end
end