NDinfFile = dir('out*.txt');
% export recorded data from te ND moview into a data file (txt) 
% copy file contents and paste into another text file
% delete the first one

lineCount = 1;
nEv = 1;

for i = 1:numel(NDinfFile)
    
    NDdata = fopen(NDinfFile(i).name);
    a=fgetl(NDdata);
    while (1)
        a=fgetl(NDdata);
        if a == -1, break; end;
        userPos = strfind(a,'User');
        if ~isempty(userPos)
            eventText(nEv,1) = {lineCount};
            eventText(nEv,2) = {a(userPos(1)+9:userPos(2)-1)};
            nEv = nEv + 1;
        else
            zData2 = sscanf(a,'%f',[1,6]);
            zCoor(lineCount)= zData2(4);
            PFSoffset(lineCount)= zData2(6);

            %pause
            lineCount = lineCount+1;
        end
    end
    fclose(NDdata);
end
acqTime = '23ms';
PFSoffset = PFSoffset/1000;
PFSoffset = PFSoffset-mean(PFSoffset);
zCoor = zCoor - mean(zCoor);
plot(zCoor,'r')
hold on
zCorrected = zCoor+PFSoffset;
plot(zCorrected,'b')
plot(PFSoffset,'g')
hold off
save('zData','zCoor','zCorrected','PFSoffset')
mkdir stats
zRange = max(zCoor)-min(zCoor);
set(gcf,'Color',[1 1 1]); title(sprintf('z position. z displacement range:%.02fum. acqTime:%s',zRange,acqTime)); xlabel('frames'); ylabel('Z');
legend('zCoor','zCorrected','PFSoffset')
imgFig = getframe(gcf); imgOut = imgFig.cdata;
imwrite(imgOut,'stats\Zcorrection.tif')