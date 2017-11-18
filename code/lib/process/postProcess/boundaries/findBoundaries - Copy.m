function [Bt,Lt,rb] = findBoundaries(R,mnFilt,t)
    rb=[];
    
    mnSz=mnFilt(1);
    mnDensity=0;
    if numel(mnFilt)>1
        mnDensity=mnFilt(2);
    end
    %% FIND BOUNDARIES
    cnvSz = 5;
    convWin = ones(cnvSz);
    Rc = conv2(R,convWin,'same'); % convolved R
    [B,L]=bwboundaries(Rc,4,'noholes');
    
    % remove small structures
    BB=[];
    for i = 1:numel(B)
        BB = [BB ; B{i}];
    end
    bimg = zeros(size(R));
	bimg(sub2ind(size(R),BB(:,1),BB(:,2)))=1;
    bimgCv = conv2(bimg,ones(3),'same');
    L(bimgCv>0)=0;
    %imagesc(L);
    szSt=hist(L(:),max(L(:))+1);
    szSt(1)=[];
    ixDel1 = find(szSt<mnSz);
    
    % remove sparse rec. structs
    LR = L.*im2bw(R,0);
    nrSt=hist(LR(:),max(L(:))+1);
    nrSt(1)=[];
    densSt=nrSt./szSt;
    ixDel2 = find( densSt<mnDensity );
    ixDel = [ixDel1 ixDel2];
    B(ixDel) = [];
    L(ismember(L,ixDel))=0;
    

    %% FIND TIGHT BOUNDARIES    
    dr=[];Bt={};s=1;Lt=[];
    %if numel(B)>1, disp('multiple sections'); end;
    %hw = waitbar(0,'sectioning');
    R0 = R;
    ixSt = unique(L);
    ixSt(1)=[];
    for k = 1:length(B) % all sections
        %if k~=19,continue; end;
        boundary = B{k};
        R = R0;
        R(L~=ixSt(k))=0;
    
        %%
        cvSz=3;
        Rinv = im2bw(R,0);
        Rinv = 1-Rinv;
        Rcv = conv2(Rinv,ones(cvSz),'same');
        Rcv(Rcv<cvSz^2)=0;
        Rcv(Rcv>0)=1;
        Rcv(1,:)=1;Rcv(end,:)=1;Rcv(:,1)=1;Rcv(:,end)=1;


        % remove internal space
        [~,Lt00]=bwboundaries(Rcv,'noholes'); % tight boundaries
        Rcv2 = Rcv;
        Rcv2(Lt00>1)=0;
        Rcv3 = conv2(Rcv2,ones(cvSz),'same');
        Rcv3(Rcv3>0)=1;
        if rem(cvSz,2)==0, Rcv3=circshift(Rcv3,[1 1]); end
        Rcv3 = 1-Rcv3;

        [Bt0,Lt0]=bwboundaries(Rcv3,'noholes'); % tight boundaries

        for i = 1:numel(Bt0) %find each tight section area
            Lt0_ = zeros(size(Lt0));
            Lt0_(Lt0==i)=1;
            if sum(Lt0_(:)) < mnSz, continue; end
            Lt(:,:,s) = Lt0_;
            Bt = [Bt;Bt0(i)];
            s = s + 1;
        end
        if 0
            figure(1); imagesc(R); axis image;
            figure(2); imagesc(Rinv); axis image;
            figure(3); imagesc(Rcv); axis image;
            figure(4); imagesc(Lt00); axis image;
            figure(5); imagesc(Rcv2); axis image;
            figure(6); imagesc(Rcv3); axis image;
            %saveFigImgFN = 'findBoundariesDBG.tif';saveFigImg
            figure(19); plotBoundaries(Bt,sum(Lt,3).*R,1); axis image;
        end
    end
    
end






