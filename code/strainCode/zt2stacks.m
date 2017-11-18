% reads fname4D (XYZT) and write separate XYZ files

nZ = StackNum; % # of z planes
imageInfo=imfinfo(fname4D);
numFrames=length(imageInfo);
filenameMaxProj_out = 'maxProj.tif';
filenameStack_out = 'stack_';

for i = 1 : numFrames
    A(:,:,i) = imread(fname4D,i);
end
isCrop = 0;
yy = 1:250;
xx = 1:250;
h = waitbar(0,'slicing data stack...');
for i = 1 : size(A,3)/nZ
    for j = 1 : nZ
        if isCrop, temp = A(yy,xx,(i-1)*nZ+j); else, temp = A(1:end,1:end,(i-1)*nZ+j); end;
        %Zstack(:,:,j) = temp(200:400,150:350);
        Zstack(:,:,j) = temp;
    end
    fn = sprintf('%s%02i.tif', filenameStack_out,i);
    stackWrite(Zstack,fn); 
    maxProj(:,:,i) = max(Zstack,[],3);
    waitbar(i / size(A,1)*nZ)
end
stackWrite(maxProj,filenameMaxProj_out);
fname = filenameMaxProj_out;
save('fname','fname');
close(h)
