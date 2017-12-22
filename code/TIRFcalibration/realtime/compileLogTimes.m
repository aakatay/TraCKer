clear all;
load cfgRT
fclose('all');
cfg = cfgSave;
waWin = cfg.waWin;
stdWin = cfg.stdWin;

fid1 = fopen(cfg.logWA);
fid2 = fopen(cfg.logThresh);
fid3 = fopen(cfg.logPos);
fid4 = fopen(cfg.logTrace);
fid5 = fopen(cfg.logSNR);


dn = [0 waWin-1 waWin-1 waWin-1 waWin ]; % frame shifts
for i = 1:5
    eval(sprintf('fid=fid%i;',i))
    fline = fgetl(fid);
    [ms] = sscanf(fline,'start time m= %i secs=%f');
    m = ms(1); s = ms(2);
    MS(i) = m*60+s;
    dMS = MS(i)-MS(1);
    while ~feof(fid)
        txt = fread(fid,11);
        if strcmp(char(txt'),'updated    ')
            fline = fgetl(fid);
            [nt] = sscanf(fline,'n=%i time=%f');
            n = nt(1); t = nt(2);
            loopT(n+dn(i),i) = t+dMS;
        else % dump the line
            fgetl(fid); 
        end
    end
end
save('timeSync','loopT')