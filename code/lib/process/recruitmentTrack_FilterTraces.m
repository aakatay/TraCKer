%% delete traces (brightness and trace length and spread )
nt = size(trInf,1); % number of traces
nd = numel(xx); % data points
delSpreadTrace = find(trInf(:,7) > thrSpread); % number of traces tb deleted
ndelSpreadTrace = numel(delSpreadTrace);
disp(sprintf('removing spread traces: %.02f%% of the traces are deleted',(ndelSpreadTrace)/size(trInf,1)*100 ));

i=1;
hw = waitbar(0);
while i <= size(trInf,1)
    if ~isempty(find(i==delSpreadTrace))
        ixdel = [trInf(i,3)-1 trInf(i,3)+trInf(i,2)];
        trInf(i+1:end,3) = trInf(i+1:end,3) - trInf(i,2);
        if ixdel(1) == 0
            xx = xx(ixdel(2):end);
            yy = yy(ixdel(2):end);
            int = int(ixdel(2):end);
            fr = fr(ixdel(2):end);
            trInf = trInf(i+1:end,:);
        elseif ixdel(2) == numel(xx)
            xx = xx(1:ixdel(1));
            yy = yy(1:ixdel(1));
            int = int(1:ixdel(1));
            fr = fr(1:ixdel(1));
            trInf = trInf(1:i-1,:);
        else
            xx = [xx(1:ixdel(1)) xx(ixdel(2):end)];
            yy = [yy(1:ixdel(1)) yy(ixdel(2):end)];
            int = [int(1:ixdel(1)) int(ixdel(2):end)];
            fr = [fr(1:ixdel(1)) fr(ixdel(2):end)];
            trInf = [trInf(1:i-1,:); trInf(i+1:end,:)];
        end
        delSpreadTrace = delSpreadTrace-1;
    else
        i = i + 1;
    end
    waitbar(i/size(trInf,1),hw,'recruitmentTrack: filtering out traces...')
end
close(hw)
disp(sprintf('removing traces: %.02f%% of the data points deleted',(nd-numel(xx))/nd*100 ));