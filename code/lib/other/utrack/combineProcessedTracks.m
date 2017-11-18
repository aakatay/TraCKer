% make section numbers to have same number of digits by adding preceding
% zeros
fn = rdir('*\**\ProcessedTracks.mat');
numFN = numel(fn); % number of sections
numFN = 30;

sectionlength = 200;
tracksAll = [];
hWB = waitbar(0,'processing sections ...')
for i = 1:numFN
    load(fn(1).name);
    for j = 1:size(tracks,2)
        tracks(j).f = tracks(j).f + (i-1)*sectionlength; 
    end
    tracksAll = [tracksAll tracks];    
    waitbar(i/numFN,hWB,'processing sections ...');
end
close(hWB);

tracks = tracksAll;
save ProcessedTracks.mat tracks -v7.3