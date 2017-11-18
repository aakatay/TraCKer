function [numFrames] = findNumFrames(fname)

    err = 1;
    load fname

    N=19;
    n2=N+1;
    frm2 = 2^n2;
    frm1=1;
    found = 0;
    numFrames =0;
    while ~found 
        while ~isempty(err)
            try 
                err = [];
                [~] = imread(fname,uint16(frm1+frm2/2));
            catch err
                frm2 = round(frm2/2);
            end
        end
           
        frm1 = frm1+frm2/2;
        %frm2 = frm1+frm2;
        err = 1;
        if frm2 == 2
            found = 1;
            numFrames = frm1;
        end
    end 
end