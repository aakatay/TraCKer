%fCoeff = dir('Coeff*');
if exist(coeffMat), load(coeffMat); end;
load fname; 
nCh = numel(fname);
for i = 1: nCh
    if fname(i) == '_'
        last_ = i;
    end
end

if ~exist('last_')
    label = [];
else
    label= fname(last_:end-4);
end
traceJmplessDataFileNm = sprintf('traceJmplessData-coeff%d%s.mat',round(Coeff),label);
traceDataFileNm = sprintf('traceData-coeff%d%s.mat',round(Coeff),label);
traceDataFileNm0 = sprintf('traceData0-coeff%d%s.mat',round(Coeff),label);
posDataFileNm = sprintf('posData-coeff%d%s.mat',round(Coeff),label);
xyzDataFileNm = sprintf('xyzData-coeff%d%s.mat',round(Coeff),label);
xyzDataGausFileNm = sprintf('xyzDataGaus-coeff%d%s.mat',round(Coeff),label);
xyzDataGausFiltFileNm = sprintf('xyzDataGausFilt-coeff%d%s.mat',round(Coeff),label);
spotWinFileNm = sprintf('spotWin-coeff%d%s.mat',round(Coeff),label);
nbinsFileNm = sprintf('nbins-coeff%d%s.mat',round(Coeff),label);
tracePlotFileNm = sprintf('plotData-coeff%d.mat',round(Coeff));
%cfgLocalizeFileNm = sprintf('cfgLocalize-coeff%d.mat',round(Coeff));
cfgTraceFileNm = sprintf('cfgTrace-coeff%d%s.mat',round(Coeff),label);
cfgGausFitFileNm = sprintf('cfgGausFit-coeff%d%s.mat',round(Coeff),label);
nmCoeff = sprintf('Coeff%s.mat',label);