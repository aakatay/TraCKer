
if exist('BREAK.MAT' ,'file')
    if exist('breakVal.mat','file')
        clear
        load breakVal;
        delete('BREAK.MAT','breakVal')
    else
        save('breakVal');
        return;
    end
end
