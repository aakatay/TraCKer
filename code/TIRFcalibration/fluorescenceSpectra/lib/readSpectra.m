function s = readSpectra(fn)
    fid = fopen(fn);
    s = fscanf(fid,'%f');   
end