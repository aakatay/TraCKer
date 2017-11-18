% converts the CME results array indexing to linear indexing for
% recruitmentTrack analysis
% BEFORE: copy CFG.mat to main directory
% AFTER: 
% - generate fname.mat with fname= ' ... .tif'
% - move output .mat files and fname.mat to '_' folder and run recruitmentTrack 
clear;


omd = 'orig_movies';
tmpd = dir(fullfile(omd,'*.tif'));
orig_movies = cell(length(tmpd),1);
for i = 1:length(tmpd)
    orig_movies{i} = fullfile(omd,tmpd(i).name);
end

smd = 'split_movies';
tmpd = dir(smd);
tmpd = tmpd([tmpd.isdir]);
tmpd(strncmp({tmpd.name},'.',1)) = [];
movies = length(tmpd);
dirname = cell(movies,1);
for i = 1:movies
    dn = fullfile(smd,tmpd(i).name);
    dirname{i} = dn;
    sectDir = rdir([dn '\section*']);
    sections(i) = numel(sectDir);
end
clear sectDir dn smd;
%create paths to all the data
SectionSize = zeros(movies,1);
moviefol = cell(max(sections),movies);
moviename = cell(max(sections),movies);
paths = cell(max(sections),movies);
for i = 1:movies
    for i2 = 1:sections(i)
        tmpn = fullfile(dirname{i},['Section',num2str(i2)]);
        tmpd = dir(tmpn);
        moviefol{i2,i} = fullfile(tmpn,tmpd(3).name,'ch1');
        paths{i2,i} = fullfile(moviefol{i2,i},'Tracking','ProcessedTracks.mat');
        if ~exist(paths{i2,i},'file')
            tracks = [];
            processingInfo = [];
            save(paths{i2,i},'tracks','processingInfo')
        end
        tmpd = dir(fullfile(moviefol{i2,i},'*.tif'));
        moviename{i2,i} = fullfile(moviefol{i2,i},tmpd.name);
        if i2==1
            SectionSize(i) = length(imfinfo(moviename{i2,i}));
        end
    end
end
%%

TraceX = [];
TraceY = [];
frmNoTrace = [];
TraceINT = [];
TraceBG = []; % background
lastFrame = 0;
NT = 0; 
trInf(1,3) = 1;
lc = 0; % cumulative l traces
for i = 1:movies % for each orig movie
    for j = 1:sections(i)
        lastFrame = SectionSize*(j-1);
        load(cell2mat(paths(j,i)),'tracks');
        for k = 1:numel(tracks) % for each trace
            ff = tracks(k).f;
            xx = tracks(k).x;
            yy = tracks(k).y;
            AA = tracks(k).A;
            cc = tracks(k).c;
            tsc = find(isnan(ff)); % trace split coordinates
            te = [tsc-1 numel(ff)]; % trace ends
            for l = 1:numel(te) % for each trace including crossing traces
                if l==1, ts=1; else ts=te(l-1)+2; end
                f = ff(ts:te(l)); 
                x = xx(ts:te(l));
                y = yy(ts:te(l));
                A = AA(ts:te(l));
                c = cc(ts:te(l));
                TraceX = [TraceX x];
                TraceY = [TraceY y];
                frmNoTrace = [frmNoTrace f+lastFrame];
                TraceINT = [TraceINT A];
                TraceBG = [TraceBG c];
                tl = numel(x); % trace length
                trInf(k+l+lc-1,1) = f(1)+lastFrame;
                trInf(k+l+lc-1,2) = tl;
                trInf(k+l+lc-1+1,3) = trInf(k+l+lc-1,3)+tl;
                if trInf(k+l+lc-1+1,3)==0
                    cccc=1;
                end
                trInf(k+l+lc-1,6) = sum(A(A>0))/tl; % intensity average
                
                fff = frmNoTrace(trInf(k+l+lc-1,3):trInf(k+l+lc-1,3)+tl-1);
                if ~isempty(find(fff(2:end)-fff(1:end-1)<0))
                    ccc=1;
                end
            end
            lc = numel(te)-1 + lc;
        end
        lc = size(trInf,1)-1;
    end
    trInf = trInf(1:end-1,:);
    jumpEl = find(TraceINT<0);
    TraceX(jumpEl) = 0;
    TraceY(jumpEl) = 0;
    TraceINT(jumpEl) = 0;
    
    
    ixSptFrm = zeros(max(frmNoTrace)+1,1); 

    fnXYZtraceData0 = 'traceData0-coeff-CMEanaly.mat';
    fnXYZ = 'xyzDataGaus-coeff-CMEanaly';
    %load CFG; % for minXYspread
    save(fnXYZtraceData0,'frmNoTrace','TraceX','TraceY','TraceINT','trInf');
    save(fnXYZ,'ixSptFrm');
    %trInf
            % 1: 1st frame
            % 2: number of frames
            % 3: position in the trace array
            % 4-6: mean x, y , int
            % 7 : std deviation from the center
    mkdir('_')        
    fname = sscanf(cell2mat(orig_movies(1)),'orig_movies\\%s');
    save('fname','fname');

    movefile('fname.mat','_\fname.mat');
    movefile('xyzDataGaus-coeff-CMEanaly.mat','_\xyzDataGaus-coeff-CMEanaly.mat');
    movefile('traceData0-coeff-CMEanaly.mat','_\traceData0-coeff-CMEanaly.mat');
    recruitmentTrack

end