function stackWriteColor(img3D,fileName)
% writes tiff stack
    sizeImg3D = size(img3D);
    imwrite(uint16(img3D(:,:,:,1)),char(fileName),'tiff','Compression', 'none');
    
    h = waitbar(0,'Writing to disk...');
    for z = 2:size(img3D,4)
        try
            imwrite(uint16(img3D(:,:,:,z)),char(fileName),'tiff','WriteMode','append','Compression', 'none');
        catch
            q=0;
            while q == 0
                try 
                    imwrite(uint16(img3D(:,:,:,z)),char(fileName),'tiff','WriteMode','append','Compression', 'none');
                    q = 1;
                catch
                    q=0;
                    disp(z)
                end
            end
        end
        waitbar(z/size(img3D,3));
    end
    close(h);
end