% reads structure data outputs from:
% 1-findEdgeDistanceDistribution.m
% 2-recruitComp.m


%% 1: thickness and sharpness(edge recr.)
edgeFNdir = rdir('edgeDistData*.mat');
[edgeFNdir,ixs1]=sortFilenameIx(edgeFNdir); % sort acc. to indices
st = struct;
for i = 1:numel(edgeFNdir)
    edgeFN = edgeFNdir(i).name;
    load(edgeFN); % matARR : [xcXLS' NdXLS' xcr' Nr' xcraXLS' NraXLS']
    A = matARR; clear matARR;
    st(i).thickness = A(end,3);
    minSz = 3;
    nrec = A(:,2);
    nrec(isnan(nrec))=0;
    st(i).nrecs = sum(nrec); % # of recs
    st(i).stix = ixs1(i);
    if size(A,1)<=minSz % small
        st(i).sharpness = nan;
        st(i).stChange = nan;
        continue
    end
    A = A(1:end,:); % remove <2px (if convWinSz:5)
    A(isnan(A))=0;
    st(i).sharpness = (sum(A(:,3).*A(:,4))-sum(A(:,1).*A(:,2)))/sum(A(:,2));
    st(i).stChange = (sum(A(:,5).*A(:,6))-sum(A(:,1).*A(:,2)))/sum(A(:,2)); % structure change
end


%% 2: recruitment position in intensity 
intFNdir = rdir('data*.mat');
%if strcmp(intFNdir(end).name,'intPlotData.xls'), intFNdir = intFNdir(1:end-1); end;
[intFNdir,ixs2]=sortFilenameIx(intFNdir); % sort acc. to indices
findMissingIx(ixs1,ixs2); % check indices of files & find any mismatch

for i = 1:numel(intFNdir)
    intFN = intFNdir(i).name;
    load(intFN); % 'ixpx','intRec','intST','th'
    A = [ixpx intRec intST];
    rec = A(1:end,2);
    int = A(1:end,3);
    recCV = smooth(rec,5);
    [~,ix]=max(recCV);
    st(i).intMxRec = int(ix);
end

%% intensity difference distribution (recruitmentImage - preImage  )
rintFNdir = rdir('rint*.mat');
[intFNdir,ixs3]=sortFilenameIx(intFNdir); % sort acc. to indices
findMissingIx(ixs1,ixs3); % check indices of files & find any mismatch

for i = 1:numel(rintFNdir)
    rintFN = rintFNdir(i).name;
    load(rintFN); % matARR : [xc' N']
    A = matARR; clear matARR;
    nrec = A(1:end,1);
    npx = A(1:end,2);
    st(i).diffInt = sum(abs(nrec.*npx))/sum(npx);
end

st = st';
aa=struct2dataset(st);
%% compile
thickness=double(aa(:,1));
npx=double(aa(:,2));
sharpness=double(aa(:,4));
stChange=double(aa(:,5));
intMxRec=double(aa(:,6));
diffInt=double(aa(:,7));
statFN = 'strStats.xls';
xlswrite(statFN,{'_#','npx', 'thickness', 'sharpness', 'intMxRec', 'diffInt', 'stChange'})
xlswrite(statFN,[ixs2' npx thickness sharpness intMxRec diffInt stChange],1,'A2')

%%
