% pulls out the recruitment stats 
clear all;
close all;
F = findall(0,'type','figure','tag','TMWWaitbar'); delete(F)

isDyn = 1;

tf = 5; % time resolution  for joiing data with diff acqTime [ms]
 % time convolve windonw to find dynamin peak [ms]
tProf = [40 20]; % dynamin profile offset [sec]
minPk = 6; %  peak intensity (after convolution) for pit alignment  

hrIT = 20; % high resolution histogram time boundary 
hrIT2 = 1/2;
hrDur = 50;
hrDur2 = 10;
if isDyn
    %mxTime = 50;
    cellType = 'DYN';
    tConv = 500; % convolution window for profile
else
    %mxTime = 120;
    cellType = 'AP2';
    tConv = 2000; % convolution window for profile
end
tc = tConv / tf; % # of frames in the convolve win

pitCoorFN = rdir('**\pitCoors.mat'); % XY LTpitsto the pits (CMEpitFinder.m)
traceFN = rdir('**\traceData_recTrack.mat'); % recruitmenTrack.m
if numel(traceFN)~=numel(pitCoorFN), error('cell data missing'); end

% output filename
CMEpitRecruitmentsFN = 'CMEpitRecruitments'; % MAT
recruitmentResultsFN = 'recruitmentResults'; % TIF
recruitmentTraceResultsFN = 'recruitmentTraceResults'; % TIF
    

%% pit info
nCells = numel(traceFN);
trInf2 = [];
frstFrame = [];
lastFrm = [];
acqT = [];
cellNo = [];
duration = [];
NP = 0;
J = 0;
RDATA = nan(1,40000);
RINTdisp = nan(1,40000);
for i = 1:nCells
    load(pitCoorFN(i).name);
    load(traceFN(i).name)

    np = size(XY,1); % number of pits
    XY = round(XY/4);
    x = XY(:,1);
    y = XY(:,2);

    rX = trInf(:,4); % recruitment data
    rY = trInf(:,5);

    rData =[];
    rintdisp = [];
    % output file
    for j = 1:np % each pit
        Lt = LtPits(:,:,j);
        LtWide = conv2(Lt,ones(3),'same');
        ixs = find(Lt(sub2ind(size(Lt),round(rY*4),round(rX*4))));
        ixs = find(LtWide(sub2ind(size(Lt),round(rY*4),round(rX*4))));
        fr = trInf(ixs,1);
        
        % recruitment stats
        nr = numel(ixs);
        nfr = trInf(ixs,2); % number of frames
        intRec = trInf(ixs,11); % number of frames
        for k = 1:nr
            rData(j,fr(k),1) = intRec(k); % 
            rData(j,fr(k),2) = nfr(k); % 
            rintdisp(j,fr(k):fr(k)+nfr(k)) = intRec(k); % 
        end
        
        intTime = fr - circshift(fr,1);
        intTime(1)=[];
        ITfrm(J+j,1:numel(intTime)) = intTime; % interval time
        IT(J+j,1:numel(intTime)) = intTime*acqTime*1e-3; % interval time
        FR(J+j,round(fr*acqTime/tf)) = 1; % recruitment times (binary)
        
        trInf2_ = trInf(ixs,:);
        trInf2 = [trInf2; trInf2_];
        %strInt = sprintf('%s,''No%i,int:%i''',strInt,i,Int(i));
        acqT =  [acqT acqTime];
        cellNo = [cellNo i];
        % duration
        frstFrame = [frstFrame min(fr)];
        lastFrm = [lastFrm max(fr)];
        duration =  [duration acqTime*(lastFrm(end)-frstFrame(end)+1)*1e-3]; % [sec]
        
    end
    RDATA(NP+1:NP+np,1:size(rData,2),1:2)=rData;
    RINTdisp(NP+1:NP+np,1:size(rintdisp,2))=rintdisp;
    NP = NP + np;
    J = J + np;
end
% output filename
CMEpitRecruitmentsFN = [CMEpitRecruitmentsFN '-' cellType '_' num2str(NP) 'pits_' num2str(nCells)  'cells.mat'];
recruitmentResultsFN = [recruitmentResultsFN '-' cellType '_' num2str(NP) 'pits_' num2str(nCells)  'cells.tif'];
recruitmentTraceResultsFN = [recruitmentTraceResultsFN '-' cellType '_' num2str(NP) 'pits_' num2str(nCells)  'cells.tif'];
    
save(CMEpitRecruitmentsFN,'frstFrame','lastFrm','duration','FR','cellNo','acqT')

%% =========== rint TIF STACK =========================
delete(recruitmentTraceResultsFN);
%% recs intensity
figure(5000)
rint = RDATA(:,:,1); % intensity
rfrm = RDATA(:,:,2); % number of frames

lastRec = find(sum(~isnan(rint),1),1,'last');
rint = rint(:,1:lastRec);
rint(isnan(rint))=0;

rintDisp = RINTdisp(:,1:lastRec);
rintDisp(isnan(rintDisp))=0;
imagesc(rintDisp)

% intensity hist
figure(5001)
mxInt = 4e4;
rintHist = rint(:);
rintHist(rintHist==0)=[];
ixIntFault = find(rintHist>mxInt);
ixIntFaultTx = sprintf('%.01f%% of recs extreme intensity(>%i) and removed',100*numel(ixIntFault)/numel(rintHist),mxInt);
rintHist2 = rintHist;
rintHist2(ixIntFault) = [];
hist(rintHist2,100);
ylabel('counts')
xlabel('intensity')

    set(gcf,'color','w');  
    set(gca, 'box','off')
save('rintHist2','rintHist2')
grid minor
title(ixIntFaultTx)
        imgFig = getframe(gcf); 
        figCap = imgFig.cdata;
        imwrite(figCap,recruitmentTraceResultsFN,'WriteMode','append','Compression', 'none') 


%% =========== PITS TIF STACK =========================
delete(recruitmentResultsFN);
%% display pit times
figure(501)
t = 1:size(FR,2)*tf*1e-3;
npits = 1:size(FR,1);
FRc = conv2(FR,ones(1,tc));
imagesc(t,npits,FRc)
ylabel('pits');
xlabel('seconds')
title(sprintf('recruitments %ims convolution',tConv))
        imgFig = getframe(gcf); 
        figCap = imgFig.cdata;
        imwrite(figCap,recruitmentResultsFN,'WriteMode','append','Compression', 'none') 

%% intervl times
figure(502)
subplot(2,2,1)
IT(IT==0) = nan;
hist(IT(:),100)
xlabel('seconds'); title('interval histogram')

subplot(2,2,2)
IT2 = IT;
IT2(IT2>hrIT) = nan;
hist(IT2(:),100)
xlabel('seconds'); title('interval histogram')

subplot(2,2,3)
IT3 = IT;
IT3(IT3>hrIT2) = nan;
hist(IT3(:),20)
xlabel('seconds'); title('interval histogram')

subplot(2,2,4)
nf=20;
ITfrm2 = ITfrm;
ITfrm2(ITfrm2==0) = nan;
ITfrm2(ITfrm2>nf)=nan;
hist(ITfrm2(:),nf-1)
xlabel('frames'); title('interval histogram')

        imgFig = getframe(gcf); 
        figCap = imgFig.cdata;
        imwrite(figCap,recruitmentResultsFN,'WriteMode','append','Compression', 'none') 

%% pit duration
figure(503)

subplot(2,2,1)
duration2 = duration;
duration2(duration2>hrDur) = [];
hist(duration2,20)
xlabel('seconds'); title('recruitment duration histogram')

subplot(2,2,2)
duration3 = duration;
duration3(duration3>15) = [];
hist(duration3,20)
xlabel('seconds'); title('recruitment duration histogram')

        imgFig = getframe(gcf); 
        figCap = imgFig.cdata;
        imwrite(figCap,recruitmentResultsFN,'WriteMode','append','Compression', 'none') 

%% profile
figure(511)
[pk,ix]=max(FRc,[],2); % peak value and time
%figure; hist(pk)
ix2 = ix(find(pk>=minPk));
%ix2 = ix;
td1 = tProf(1)/tf*1e3;
td2 = tProf(2)/tf*1e3;
clear FRp;
for i = 1:numel(ix2) % each dynamin peak
    frm0 = ix(i) - td1 + td1; 
    frm2 = ix(i) + td2 + td1;
    FRc = padarray(FR,[0 td1]);
    FRp(i,:) = FRc(i,frm0:frm2);
end
FRpc = conv(sum(FRp,1),ones(1,tc),'same');

td = td1 + td2;
tt = [1:td+1]*tf*1e-3;
plot(tt,FRpc)
xlabel('seconds')
ylabel('au')
title(sprintf('recruitment profile %ims convolution (aligned at peaks)',tConv))

        imgFig = getframe(gcf); 
        figCap = imgFig.cdata;
        imwrite(figCap,recruitmentResultsFN,'WriteMode','append','Compression', 'none') 




%%
figure(515)
intP = sum(FR,2);
[int, ix]  = sort(intP,'descend');

tx = strtrim(cellstr(num2str(ix))'); % string(ix)
plot(int);
set(gca,'XTick', 1:NP);
set(gca,'XTickLabel', tx)
set(gca,'XTickLabelRotation', 90)
 
xlim([-1 NP+1])
ylim([-1 max(intP)+1])
xlabel('pit index')

ylabel('# of recruitments')
%maximize
grid 

        imgFig = getframe(gcf); 
        figCap = imgFig.cdata;
        imwrite(figCap,recruitmentResultsFN,'WriteMode','append','Compression', 'none') 

%% XLS write
xlsFile = 'CMEpitRecruitments';
xlsFile = [xlsFile '-' cellType '_' num2str(NP) 'pits_' num2str(nCells)  'cells.xls'];
delete(xlsFile)
xlswrite(xlsFile,[{'sheet2 : intervals'}; {'sheet3 : durations'}; {'sheet4 : intensity'}],1)
xlswrite(xlsFile,IT(:),2)
xlswrite(xlsFile,duration',3)
xlswrite(xlsFile,intP,4)

