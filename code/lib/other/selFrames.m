CD = cd;
cropFolders = dir('crop*');
while isempty(cropFolders) % go to the main folder
    cd ..
    cropFolders = dir('crop*');
end

% find background info
bckDataFiles = rdir('*\**\CoeffNormMulSmth.mat');

CoeffNormAll = [];
lastVal =1;
for i=1:numel(bckDataFiles) % join coeff values
    load(bckDataFiles(i).name);
    if ~isempty(CoeffNormAll)
        lastVal = CoeffNormAll(end);
    end
    CoeffNormAll = [CoeffNormAll;CoeffNormMulSmth*lastVal];
end
plot(CoeffNormAll); title('Coefficient change'); ylabel('Coefficient'); xlabel('frames');
cd(CD)