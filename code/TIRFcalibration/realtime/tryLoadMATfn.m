function out=tryLoadMATfn(evalTxt,Tloop)
% tries in loop to load MATrtDetectThresh.mat MATrtTraCKerPos.mat ...
    out=[];
    isErr = 1;
    while isErr
        isErr = 0;
        try 
            eval(evalTxt)
            if ~isfield(out,'lmpState')
                isErr = 1;
            elseif ~isfield(out,'nFrst')
                isErr = 1;
            elseif ~isfield(out,'nLast')
                isErr = 1;
            end
        catch 
            isErr = 1;
            pause(Tloop);
        end
    end
end
