% extracts the high layers from the data to determine the outline of the
% cell
clear;
fname = rdir('*.nd2')


for i = 1:size(fname,1)
    fn(i) = str2num(strcat(fname(i).name(5),fname(i).name(1:3)));
end
[Y,J] = sort([fname.datenum]); 

for i = 1:size(fname,1)
    ix = J(i);
    file_out = strcat(sprintf('%02d',i),'-',fname(ix).name(1:3),'.tif');
    data = bfopen(cell2mat(fn(i))cell2mat(fn(i)));
    metadata = data{1, 2};
    A = data{1,1};
    for j = 1 : size(A,1)
        maxProj(:,:,j) = cell2mat(A(j,1));
    end
    stackWrite(maxProj,file_out);
end
