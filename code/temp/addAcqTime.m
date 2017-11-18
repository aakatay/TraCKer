% add acqtime
fid = fopen('frames.txt');
gl = fgetl(fid);
gl = fgetl(fid);
acqTime = str2num(gl(9:end));
fclose(fid);

save('traceData_recTrack','acqTime','-append')