% browse to the root folder
ixRt = 1;
while 7~=exist('code') || 7~=exist('data') % not in the main folder
    cd('..');
    ixRt = ixRt+1;
    if ixRt>10
        display('manually change the current directory to DSrotate main folder');
        break;
    end
end 
if ~isempty([strfind(cd,'code') strfind(cd,'data')])
    cd('..')
end