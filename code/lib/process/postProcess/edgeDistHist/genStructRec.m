% generates structures to compare the results from
% findEdgeDistanceDistribution.m

% REF 1 : uniform distribution
[sy,sx] = size(R); % size
sg = ones(sy,sx); % flat structure
sg = sg.*Lt(:,:,k); % struct generated
%figure(1);imagesc(sg)
[Y,X,N]=findPos(sg);

stSz = numel(Y); % structure size
drmin=[];
j1=1;
for i2 = 1:numel(Y)
   j2 = j1 + N(i2)-1;
   dr_ = sqrt((boundaryT(:,1)-Y(i2)).^2+(boundaryT(:,2)-X(i2)).^2);
   drmin(j1:j2) = min(dr_);
   j1 = j2+1;
end

if isdrAmin
    % REF 2 : pre bleach image
    sg = repelem(double(imread(fnameA)),4,4); % pre bleach image

        %remove background
        BACKmean=[min(mean(sg,1)),min(mean(sg,2))];
        BACK =min(BACKmean);
        sg = sg-BACK;
        sg(sg<0)=0;

    sg = round(sg);
    sg = sg.*Lt(:,:,k); % struct generated
    [Y,X,N]=findPos(sg);
    drAmin=[];
    j1=1;
    for i3 = 1:numel(Y)
        j2 = j1 + N(i3)-1;
       dr_ = sqrt((boundaryT(:,1)-Y(i3)).^2+(boundaryT(:,2)-X(i3)).^2);
       drAmin(j1:j2) = min(dr_);
        j1 = j2+1;
    end
end

