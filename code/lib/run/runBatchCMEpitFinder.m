% batch processes cells
% run following before: 
%   Coeff: TraCKer_3D_w_ZcolorPlot_deep.m
%   prebleach: CMEpitFinder.m
%   selPitMean*.tif: CMEpitImage(isTL = 1) (requires manual cropping)
% update following:
%   acqTime in recruitmentTrack.m
%   isDyn in CMEpitFinder.m (not needed anymore)
%   isDyn in CMEpitRecruitments.m
%   isTL in CMEpitImage.m
%   cellType in compileCMEpitImage.m

% later check
% bleach half time : recruitmentRate-hlfTimeXXXsec


ix = [3, 5, 6, 7, 8, 12]; % cell indices
ix = [8, 12]; % cell indices
save('cellIX','ix')
for  i = 1:numel(ix)    
    load('cellIX')
    fn = sprintf('cell%i',ix(i));
    cd(fn)
    TraCKer_3D_w_ZcolorPlot_deep
    recruitmentTrack
    addAcqTime% temporary
    CMEpitFinder
    CMEpitImage
    cd('..')
end


CMEpitRecruitments; % statistical recruitment analysis
compileCMEpitImage