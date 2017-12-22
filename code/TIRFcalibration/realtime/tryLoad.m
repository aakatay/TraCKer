function out=tryLoad(evalTxt,Tloop)
% tried in loop
    out=[];
    isErr = 1;
    while isErr
        isErr = 0;
        try 
            eval(evalTxt)
        catch 
            isErr = 1;
            pause(Tloop);
        end
    end
end
