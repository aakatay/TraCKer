function fid = safefopen(fileName)
   fid = fopen(fileName, 'w');
   c = onCleanup(@()fclose(fid));

   %imread('adfasdfasd');
end   % c executes fclose(fid) here