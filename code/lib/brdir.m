% back recursive directory search
function folder = brdir(Fname)
% Fname : folder name

    folder=Fname;
    if isempty(dir(folder)), folder=(['../' Fname]); end;
    if isempty(dir(folder)), folder=(['../../' Fname]); end;
    if isempty(dir(folder)), folder=(['../../' Fname]); end;
    if isempty(dir(folder)), folder=(['../../../' Fname]); end;
    if isempty(dir(folder)), folder=(['../../../../' Fname]); end;
    
end