clear all;
close all;
load('traceData0-coeff0_002.mat','trInf')

load('posData-coeff0_002.mat', 'Frames','En1','Boy1')
% filter traces
%minXYspread = 0.8;
ROItracesX = (trInf(:,4)>1) .* (trInf(:,4)<=En1);
ROItracesY = (trInf(:,5)>1) .* (trInf(:,5)<=Boy1); 
ROItraces = ROItracesX .* ROItracesY;
trInf = trInf(find(ROItraces),:);


trInf2 = trInf;
trInf2(trInf(:,1)~=1,:) = [];
size(trInf,1)
size(trInf2,1)

hist(trInf2(:,2),Frames)
title('trace length (frames)')

%trInf

        % 1: 1st frame
        % 2: number of frames
        % 3: position in the trace array
        % 4-6: mean x, y , int
        % 7 : std deviation from the center