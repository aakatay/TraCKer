function WAmean
    % walking average mean

    waWin = 15; % walkiong average window length

    fn=dir('cy3_*.tif');
    fname0 = fn(1).name;

    fname = [fname0(1:end-5) '2.tif'];
    save fname fname;
    if exist(fname),delete(fname); end


    iminf = imfinfo(fname0);

    frames = numel(iminf);
    Awa = zeros(iminf(1).Height,iminf(1).Width,waWin);
    for i = 1:frames
        Awa = circshift(Awa,[0 0 1]);
        Awa(:,:,1) = imread(fname0,i);
        if i>=waWin
            A = uint16(mean(Awa,3));
            ME = 1;
            while ~isempty(ME)
                ME = [];
                try imwrite(A,fname,'WriteMode','append');
                catch ME
                    disp(sprintf('error: %i',i)); %#ok<DSPS>
                end
            end
            pause(0.1)
        end
    end
    
    
end