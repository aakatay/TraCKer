% find means displacements betw. frames and saves it 
% (data used by TraCKer_3D_w_ZcolorPlot_deep to correct for vibration)
% run this code in ventral folder

if isempty(strfind(cd, 'ventral3D'))
    display('run this code in ventral folder')
end
d = rdir('traceData-*');
[Y,I] = sort([d.datenum]); 
traceDataPath = {d(I(end)).name}; % load the most recent traceData
load traceData-coeff183.mat

%% remove NAN traces
FrstFrm = TraceX(:,1);
TraceX = TraceX(~isnan(FrstFrm),:);
TraceY = TraceY(~isnan(FrstFrm),:);
TraceZ = TraceZ(~isnan(FrstFrm),:);

Nfrm = size(TraceX,2);

x = logical(TraceX);
x2 = logical(zeros(size(x)));
x2(:,2:end) = x(:,1:end-1);
x3= logical((x2).*x); % logical array refering to px where speed can be evaluated


TraceXlast = zeros(size(TraceX));
TraceXlast(:,2:end) = TraceX(:,1:end-1);
TraceXdiff = TraceX-TraceXlast;
TraceXdiff = TraceXdiff.*(x3);

TraceYlast = zeros(size(TraceY));
TraceYlast(:,2:end) = TraceY(:,1:end-1);
TraceYdiff = TraceY-TraceYlast;
TraceYdiff = TraceYdiff.*(x3);

TraceZlast = zeros(size(TraceZ));
TraceZlast(:,2:end) = TraceZ(:,1:end-1);
TraceZdiff = TraceZ-TraceZlast;
TraceZdiff = TraceZdiff.*(x3);

TraceDiff = sqrt(double(TraceXdiff.^2 + TraceYdiff.^2 + TraceZdiff.^2));

traceSpeed = 50;
% delete non static traces
[ixTrc speedTrc] = find(TraceDiff > traceSpeed );
TraceXdiff(ixTrc,:) = 0;
TraceYdiff(ixTrc,:) = 0;
TraceZdiff(ixTrc,:) = 0;
fprintf('traceSpeed: %i, # of non-static traces: %i \n',traceSpeed,length(speedTrc));

Ntrc = size(TraceX,1);

%% calculate mean displacements
displXmean = zeros(1,Nfrm);
displYmean = displXmean;
displZmean = displXmean;
vibrCorrDataMask = ones(size(TraceX)); % mask definining the vibration corrected data

for i = 1: Nfrm   % for each frame
    [ixNz tempVal] = find(TraceXdiff(:,i) > 0);
    if isempty(ixNz)
        vibrCorrDataMask(:,i)=0;
    else
        displXmean(i) = mean(TraceXdiff(ixNz,i));
        displYmean(i) = mean(TraceYdiff(ixNz,i));
        displZmean(i) = mean(TraceZdiff(ixNz,i));
    end
end

%% save the offset values (TraCKer use these data to correct for vibration)
save('vibrationOffset.mat','displXmean','displYmean','displZmean');
copyfile('vibrationOffset.mat','../dorsal3D/vibrationOffset.mat');

%% calculates the corrected ventral values for comparison
ixnzMask = find(TraceX~=0); % nonzero mask indices
nzMask = zeros(size(TraceX));
nzMask(ixnzMask) = 1;

TraceX = TraceX - repmat(displXmean,Ntrc,1);
TraceY = TraceY - repmat(displYmean,Ntrc,1);
TraceZ = TraceZ - repmat(displZmean,Ntrc,1);

TraceX = TraceX .* nzMask;
TraceY = TraceY .* nzMask;
TraceZ = TraceZ .* nzMask;

% remove the data points that can not be corrected
TraceX = TraceX.*vibrCorrDataMask;
TraceY = TraceY.*vibrCorrDataMask;
TraceZ = TraceZ.*vibrCorrDataMask;

fileTraceData = 'traceDataVibrCorrected_preview';
save_traceData

% X = X - repmat(Xvibr,Ntrc,1);
% Y = Y - repmat(Yvibr,Ntrc,1);
% Z = Z - repmat(Zvibr,Ntrc,1);