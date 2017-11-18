function save_fname(varargin)
if nargin == 1
    fn=dir('*.tif');
    fnum = num2str(cell2mat(varargin(1)), '%03i');
    fnameFound = 0;
    for i = 1:numel(fn)
        fname_=fn(i).name;
        if ~strcmp(fname_(1:3),'MAX') && ~strcmp(fname_(1:3),'AVG')
            if strcmp(fname_(end-6:end-4),fnum)
                fname = fname_;
                fnameFound = fnameFound + 1;
            end
        end
    end
    if fnameFound == 1
        save('fname','fname');
    else
        error('multiple files or no found with the same number')
    end
end