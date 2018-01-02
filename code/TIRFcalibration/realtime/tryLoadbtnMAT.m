function out=tryLoadbtnMAT(evalTxt,Tloop)
% tries in loop to load signals\btnMAT.mat 
    out=[];
    isErr = 1;
    while isErr
        isErr = 0;
        try 
            eval(evalTxt)
            if ~isfield(out,'btnSave')
                isErr = 1;
            elseif ~isfield(out,'btnStart')
                isErr = 1;
            elseif ~isfield(out,'btnSnap')
                isErr = 1;
            elseif ~isfield(out,'btnSync')
                isErr = 1;
            elseif ~isfield(out,'btnStop')
                isErr = 1;
            end
        catch 
            isErr = 1;
            pause(Tloop);
        end
    end
end
