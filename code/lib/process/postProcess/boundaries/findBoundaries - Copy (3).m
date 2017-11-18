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
    B5 = B; L5 = L;
    
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
    L5d=L;
    

    %% FIND TIGHT BOUNDARIES    
    dr=[];Bt={};s=1;Lt=[];
    %if numel(B)>1, disp('multiple sections'); end;
    %hw = waitbar(0,'sectioning');
    R0 = R;
    ixSt = unique(L);
    ixSt(1)=[];
    kk=1;
    for k = 1:length(B) % all sections
        %if k~=19,continue; end;
        boundary = B{k};
        brdrPos = sub2ind(size(R),boundary(:,1),boundary(:,2));
        crnrPos = sub2ind(size(R),[size(R,1) 1 size(R,1) 1 ],[1 size(R,2) size(R,2) 1]);
        if sum(ismember(brdrPos,crnrPos)==1), continue; end;
        
        R = R0;
        R(L~=ixSt(kk))=0;
    
        %%
        cvSz=3;
        cvSz=3;
        Rcv = conv2(R,ones(cvSz),'same');
        %Rcv(1,:)=9;Rcv(end,:)=9;Rcv(:,1)=9;Rcv(:,end)=9;
        RcvBin = im2bw(Rcv,0);
        RcvBininv = im2bw(RcvBin,0);
        RcvBininv = 1-RcvBininv;
        


        % remove internal space
        [~,Lt00]=bwboundaries(Rcv,'noholes'); % tight boundaries
        Rcv2 = RcvBininv;
        Rcv2(Lt00>1)=0;
        Rcv3 = conv2(Rcv2,ones(cvSz),'same');
        Rcv3_ = Rcv3;
        Rcv3(Rcv3>0)=1;
        if rem(cvSz,2)==0, Rcv3=circshift(Rcv3,[1 1]); end
        Rcv4 = 1-Rcv3;

        [Bt0,Lt0]=bwboundaries(Rcv4,'noholes'); % tight boundaries

        for i = 1:numel(Bt0) %find each tight section area
            Lt0_ = zeros(size(Lt0));
            Lt0_(Lt0==i)=1;
            if sum(Lt0_(:)) < mnSz, continue; end
            Lt(:,:,s) = Lt0_;
            Bt = [Bt;Bt0(i)];
            s = s + 1;
        end
        if 1
            figure(671);imagesc(im2bw(L5,0)); axis image;
            %figure(672);imagesc(im2bw(Rc,0)); axis image;
            figure(670);imagesc(im2bw(R,0)); axis image;
            figure(674);imagesc(im2bw(bimg,0)); axis image;
            figure(675);imagesc(im2bw(bimgCv,0)); axis image;
            figure(676);imagesc(L5d); axis image;
            
            
            %figure(1000); imagesc(R); axis image;
            figure(1001); imagesc(R); axis image;
            figure(1002); imagesc(Rcv); axis image;
            figure(1003); imagesc(RcvBin); axis image;
            %figure(1004); imagesc(Lt00); axis image;
            figure(1005); imagesc(RcvBininv); axis image;
            figure(1006); imagesc(Rcv3_); axis image;
            figure(1007); imagesc(Rcv3); axis image;
            figure(1008); imagesc(Rcv4); axis image;
            %saveFigImgFN = 'findBoundariesDBG.tif';saveFigImg
            figure(1010); plotBoundaries(Bt,sum(Lt,3).*R,1); axis image;
        end
        kk=kk+1;
    end
    
end






