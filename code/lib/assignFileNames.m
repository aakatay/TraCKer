%fCoeff = dir('Coeff_*');
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
Coeff_ = CoeffFit(1);
traceJmplessDataFileNm = sprintf('traceJmplessData-coeff%d%s.mat',round(Coeff_),label);
traceDataFileNm = sprintf('traceData-coeff%d%s.mat',round(Coeff_),label);
traceDataFileNm0 = sprintf('traceData0-coeff%d%s.mat',round(Coeff_),label);
posDataFileNm = sprintf('posData-coeff%d%s.mat',round(Coeff_),label);
xyzDataFileNm = sprintf('xyzData-coeff%d%s.mat',round(Coeff_),label);
xyzDataGausFileNm = sprintf('xyzDataGaus-coeff%d%s.mat',round(Coeff_),label);
xyzDataGausFiltFileNm = sprintf('xyzDataGausFilt-coeff%d%s.mat',round(Coeff_),label);
spotWinFileNm = sprintf('spotWin-coeff%d%s.mat',round(Coeff_),label);
nbinsFileNm = sprintf('nbins-coeff%d%s.mat',round(Coeff_),label);
tracePlotFileNm = sprintf('plotData-coeff%d.mat',round(Coeff_));
%cfgLocalizeFileNm = sprintf('cfgLocalize-coeff%d.mat',round(Coeff_));
cfgTraceFileNm = sprintf('cfgTrace-coeff%d%s.mat',round(Coeff_),label);
cfgGausFitFileNm = sprintf('cfgGausFit-coeff%d%s.mat',round(Coeff_),label);
nmCoeff = sprintf('Coeff_%s.mat',label);