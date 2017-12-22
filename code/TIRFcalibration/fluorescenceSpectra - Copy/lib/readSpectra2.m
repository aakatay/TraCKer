function MMM = readSpectra2(fn)
    M = csvread(fn);
    M(end+1,:)=0;
    MM = M';
    MMM=reshape(MM(:),3,2,numel(M)/6);
    MMM(:,2,:)=[];
    MMM = squeeze(MMM)';
end