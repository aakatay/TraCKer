clear all;
close all;

isOverlay = 1;
isRT = 1;
A = zeros(50);

% output fles
posRTfn = 'dbgPosRT.tif';
posFN = 'dbgPos.tif';
posOverlayFN = 'dbgPosOverlay.tif';

if isOverlay
    stackOverlay(posRTfn,posFN,posOverlayFN,0,[]);
else

    if isRT
        FN = dir('posData*.mat');
        n=  numel(FN);
        IX = [];

        for i = 1:n
            fn = FN(i).name;
            IX  = [IX; str2num(fn(end-7:end-4))];
        end
        [~,IXs] =sort(IX);
        FN = FN(IXs);

        for i = 1:n
            fn = FN(i).name;
            IX  = [IX; str2num(fn(end-7:end-4))];
            load(fn);
            xxRT=[];yyRT=[];
            ns(i) = numel(X);
            for j = 1:numel(X)
                A(round(Y(j)),round(X(j)),i) = 1;
                xxRT = [xxRT X(j)];
                yyRT = [yyRT Y(j)];
            end
            %save dbgPosRT xxRT yyRT; return;

        end
        stackWrite(A,posRTfn)
    else
        fn = dir('xyzDataGaus-*.mat');
        load(fn.name)
        ix = ixSptFrm;
        n = numel(ixSptFrm)-1;
        for i = 1:n
            k = ix(i):ix(i+1)-1;
            x=X(k);
            y=Y(k);
            xx=[];yy=[];
            ns(i) = numel(x);
            for j = 1:numel(x)
                A(round(y(j)),round(x(j)),i) = 1;
                xx = [xx x(j)];
                yy = [yy y(j)];
            end
            %save dbgPos xx yy; return;
            %imagesc(A(:,:,i)); pause(1)
        end
        stackWrite(A,posFN)
    end
end