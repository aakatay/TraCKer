clear all;
close all;


% 1: in RT folder run  isOverlay = 0; isRT = 1;
% 2: in bulk process folder run  isOverlay = 0; isRT = 0;
% 3: copy outputs to the same folder run isOverlay = 1

isOverlay = 0;
isRT = 0;
A = zeros(50);

% output fles
posRTfn = 'dbgPosRT.tif';
posFN = 'dbgPos.tif';
posOverlayFN = 'dbgPosOverlay.tif';

if isOverlay
    stackOverlay(posRTfn,posFN,posOverlayFN,0,[]); % red and green
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

        xxRT=[];yyRT=[];
        for i = 1:n
            fn = FN(i).name;
            IX  = [IX; str2num(fn(end-7:end-4))];
            load(fn);
            ns(i) = numel(X);
            for j = 1:numel(X)
                A(round(Y(j)),round(X(j)),i) = 1;
                xxRT = [xxRT X(j)];
                yyRT = [yyRT Y(j)];
            end
            %save dbgPosRT xxRT yyRT; return;

        end
        stackWrite(A,posRTfn)
        save('posRT','xxRT','yyRT')
    else
        fn = dir('xyzDataGaus-*.mat');
        load(fn.name)
        ix = ixSptFrm;
        n = numel(ixSptFrm)-1;
        xx=[];yy=[];
        for i = 1:n
            k = ix(i):ix(i+1)-1;
            x=X(k);
            y=Y(k);
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
        save('pos','xx','yy')
    end
end

%[xx' xxRT' yy' yyRT']
%[xx'-xxRT' yy'-yyRT']