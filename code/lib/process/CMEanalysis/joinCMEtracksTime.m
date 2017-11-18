function joinCMEtracksTime(smd,omd,Threshs,sections)

    %% some file management
    tmpd = dir(fullfile(omd,'*.tif'));
    orig_movies = cell(length(tmpd),1);
    for i = 1:length(tmpd)
        orig_movies{i} = fullfile(omd,tmpd(i).name);
    end
    tmpd = dir(smd);
    tmpd = tmpd([tmpd.isdir]);
    tmpd(strncmp({tmpd.name},'.',1)) = [];
    movies = length(tmpd);
    dirname = cell(movies,1);
    for i = 1:movies
        dirname{i} = fullfile(smd,tmpd(i).name);
    end
    %create paths to all the data
    SectionSize = zeros(movies,1);
    moviefol = cell(max(sections),movies);
    moviename = cell(max(sections),movies);
    paths = cell(max(sections),movies);
    for i = 1:movies
        for i2 = 1:sections(i)
            tmpn = fullfile(dirname{i},['Section',num2str(i2)]);
            tmpd = dir(tmpn);
            moviefol{i2,i} = fullfile(tmpn,tmpd(3).name,'ch1');
            paths{i2,i} = fullfile(moviefol{i2,i},'Tracking','ProcessedTracks.mat');
            if ~exist(paths{i2,i},'file')
                tracks = [];
                processingInfo = [];
                save(paths{i2,i},'tracks','processingInfo')
            end
            tmpd = dir(fullfile(moviefol{i2,i},'*.tif'));
            moviename{i2,i} = fullfile(moviefol{i2,i},tmpd.name);
            if i2==1
                SectionSize(i) = length(imfinfo(moviename{i2,i}));
            end
        end
    end
    
    for i = 1:movies
        for j = 1:sections(i)
            load(paths(j,i),'tracks');
            tracks.f
        end
    end
    
 
end