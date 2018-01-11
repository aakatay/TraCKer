function deleteALLparallelPools
    n = 0;
    while ~isempty(gcp('nocreate'))
        delete(gcp('nocreate'))
        disp('deleting GCP')
        n = n + 1;
    end
    disp(sprintf('%i Parallel Pools deleted',n))
end