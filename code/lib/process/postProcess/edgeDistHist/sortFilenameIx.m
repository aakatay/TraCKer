function [structFNdir,ixs] = sortFilenameIx(structFNdir)
% sort filenames in a struct acc. to their index
    for i = 1:numel(structFNdir)
        strFN = structFNdir(i).name;
        ixFN(i) = str2num(strFN(end-5:end-4));
    end
    [ixs, ixSrtFN ]= sort(ixFN);
    structFNdir=structFNdir(ixSrtFN);
end