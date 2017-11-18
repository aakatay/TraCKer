% converts the CME results array indexing to linear indexing for
% recruitmentTrack analysis
% BEFORE: copy CFG.mat to main directory
% AFTER: 
% - generate fname.mat with fname= ' ... .tif'
% - move output .mat files and fname.mat to '_' folder and run recruitmentTrack 
clear;
TreshfxycFN = rdir('orig_movies\*.mat');
load(TreshfxycFN.name); % Treshfxyc
% remove empty traces
temp = squeeze(Threshfxyc(:,2,:));
temp(temp==0)=nan;
nnn = ~isnan(temp); % defined positions
TL = sum(nnn,1); % trace lengths
NT = numel(TL); % number of traces
Threshfxyc = Threshfxyc(:,:,TL>0); % delete zero traces

frames = squeeze(Threshfxyc(:,1,:));
tracex = squeeze(Threshfxyc(:,2,:));
tracey = squeeze(Threshfxyc(:,3,:));
traceInt = squeeze(Threshfxyc(:,5,:));
clear Threshfxyc;

tracex(tracex==0)=nan;
tracey(tracey==0)=nan;
frames(frames==0)=nan;
traceInt(traceInt==0)=nan;


nnn = ~isnan(tracex); % defined positions
TL = sum(nnn,1); % trace lengths
NT = numel(TL); % number of traces


TraceX = tracex(:)';
TraceY = tracey(:)';
Frames = frames(:)';
TraceINT = traceInt(:)';
TraceX = TraceX(~isnan(TraceX));
TraceY = TraceY(~isnan(TraceY));
Frames = Frames(~isnan(Frames));
TraceINT = TraceINT(~isnan(TraceINT));

ix = 0;
jlc = 0; %jump length cumulative
trInf(1,3) = 1;
hw = waitbar(0,'CME2Tracker...');
for i = 1:NT % traces
    tl = TL(i);
     % data index
    trInf(i,1) = Frames(ix+1);
    nf = Frames(ix+tl)-Frames(ix+1)+1; % number of frames
    trInf(i,2) = nf;
    trInf(i+1,3) = trInf(i,3)+trInf(i,2);
    trInf(i,6) = mean(traceInt(~isnan(traceInt(:,i)),i)); % intensity average
    frmNoTrace(ix+jlc+1:ix+jlc+nf) = Frames(ix+1) : Frames(ix+tl);
    for j = 1:tl-1 % trace elements
        ix = ix +1;
        jl = Frames(ix+1)-Frames(ix)-1; % jump length
        if jl ~= 0 
            TraceX(ix+jl+jlc+1:end+jl) = TraceX(ix+jlc+1:end);
            TraceX(ix+jlc+1:ix+jl+jlc) = nan;
            TraceY(ix+jl+jlc+1:end+jl) = TraceY(ix+jlc+1:end);
            TraceY(ix+jlc+1:ix+jl+jlc) = nan;
            TraceINT(ix+jl+jlc+1:end+jl) = TraceINT(ix+jlc+1:end);
            TraceINT(ix+jlc+1:ix+jl+jlc) = nan;
            
        end
        jlc = jl + jlc;
    end
    ix = ix +1;
    waitbar(i/NT,hw,'CME2Tracker...');
    
end
close(hw)
ixSptFrm = zeros(max(Frames)+1,1); 

fnXYZtraceData0 = 'traceData0-coeff-CMEanaly.mat';
fnXYZ = 'xyzDataGaus-coeff-CMEanaly';
load CFG; % for minXYspread
save(fnXYZtraceData0,'cfg','frmNoTrace','TraceX','TraceY','TraceINT','trInf');
save(fnXYZ,'ixSptFrm');
%trInf
        % 1: 1st frame
        % 2: number of frames
        % 3: position in the trace array
        % 4-6: mean x, y , int
        % 7 : std deviation from the center
mkdir('_')        
fname = sscanf(TreshfxycFN.name,'orig_movies\\%s');
fname = fname(1:end-4);
save('fname','fname');

movefile('fname.mat','_\fname.mat');
movefile('xyzDataGaus-coeff-CMEanaly.mat','_\xyzDataGaus-coeff-CMEanaly.mat');
movefile('traceData0-coeff-CMEanaly.mat','_\traceData0-coeff-CMEanaly.mat');
recruitmentTrack
