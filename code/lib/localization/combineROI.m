function combineROI
% concatenates crops in a grid from a movie

    is1Dgrid = 1;
    outFN='acqPitCrops.tif';
    load pitCoors; % XY
    load fname; % acq file
    load pitCropsInfo; % frame crop info: frstFrame, lastFrm, szpad2, windowSz
    nf = mxMovieLen; % number of frames
    %nf=10;
    delete(outFN);
    
    if exist(outFN)
        disp('acqPitCrops.tif is already generated')
        return;
    else
        disp('generating acqPitCrops.tif')
    end

    % Tracker processing frames
    iminf = imfinfo(fname);
    %nf = numel(iminf);
    fo = fopen('frames.txt');
    inp = fgetl(fo);
    sc = strfind(inp,':'); % pos. semicolon
    frm1 = str2num(inp(1:sc-1));
    frm2 = str2num(inp(sc+1:end));


    % concatenation grid
    XY = round(XY/4);
    np = size(XY,1);
    if is1Dgrid
        nr = 1;
        nc = np;
    else
        nr = round(sqrt(np/2)); % # of rows
        nc = ceil(np/nr);
    end
    n = nc*nr;

    % crop sizes
    szpad = windowSz/2; % pad to read close to boundaries
    sz = szpad2*2+windowSz; % szpad2: pad to display

    NR = nr*sz;
    NC = nc*sz;

    a1d = zeros(sz,nr*nc*sz);
    fc=ones(np,1); % frame counter
    i = 1;
    ok = zeros(1,np); 
    hw=waitbar(0);
    while (1) % read frames
        A=imread(fname,frm1+i-1);
        A = padarray(A,[szpad szpad]);
        a = [];
        for j = 1:np % each crop 
            f = fc(j);
            if i<frstFrame(j), continue; end;
            if i>lastFrm(j), ok(j)=1; continue; end; 
            if fc(j)>nf, ok(j)=1; continue; end; 
            xy = XY(j,:)+szpad;
            x=xy(1); y=xy(2);
            a = A(y-szpad:y+szpad-1,x-szpad:x+szpad-1);
            a = padarray(a,[szpad2 szpad2]);
            a1d(:,1+(j-1)*sz:j*sz,f) = a; 
            fc(j)=fc(j)+1;
        end
        %a1d(:,1:size(a,2)) = a;
        if sum(ok) == np
            break; 
        end
        i=i+1;
        waitbar(i/max(lastFrm),hw,'reading frames')
    end
    close(hw)
    
    % reshape grid 1D-->2D
    hw=waitbar(0);
    for i = 1:size(a1d,3) % read frames
        a = a1d(:,:,i);
        a2d = [];
        for rix = 1:nr % each grid row 
            a2d = [a2d ; a(:,1:NC)];
            a(:,1:NC) = []; 
        end
        a2d(a2d==0)=1500; % gray frames
        %A2d(:,:,i) = a2d;
        waitbar(i/size(a1d,3),hw,'writing frames')
        imwrite(uint16(a2d),outFN,'WriteMode','append','Compression', 'none')
    end
    close(hw)
end