function [Btn,Ltn] = findBoundaries(Btn,Ltn,Rc,pt)
    % pExc : position needs to be excluded
    np = size(pt,1);
    for j = 1:np % for each selected tight request
        px = pt(j,1);
        py = pt(j,2);
        rb=[];
        pix = sub2ind(size(Rc),py,px);
        pim = zeros(size(Rc));
        pim(pix)=1;

        ns = size(Ltn,3);
        ixs = [];
        for i = 1:ns
            Lt = Ltn(:,:,i);
            if ~isempty(find(Lt.*pim>0))
                ixs = i;
                Bt = Btn{i};
                break;
            end
        end
        if isempty(ixs), continue; end;

        Ls = im2bw(sum(Lt,3),0) - im2bw(sum(Rc,3),0); % space inside structure
        Ls(Ls<0)=0;
        Ws = watershed(Ls);
        CC = bwconncomp(Ls,4);
        reg = CC.PixelIdxList; % regions
        for i = 1:numel(reg)
            if ~isempty(find(pix==reg{i}))
                break
            end
        end
        reg = reg{i};
        reg_ = zeros(size(Rc));
        reg_(reg)=1;
        reg=reg_;
        rc0 = conv2(reg,ones(2),'same');
        rc0(rc0<4)=0;
        CC2 = bwconncomp(rc0,4);
        reg2 = CC2.PixelIdxList; % regions
        for i = 1:numel(reg2)
            reg3 = reg2{i};
            rc = zeros(size(Rc));
            rc(reg3)=1;
            rcDown = circshift(rc,1,1);
            rcRight = circshift(rc,1,2);
            rcDownRight = circshift(rcRight,1,1);   
            RC = double(im2bw(rc+rcRight+rcDown+rcDownRight,0));
            if isempty(find(pix==find(RC>0)))
                reg = reg - RC;

            end
        end
        CC3 = bwconncomp(reg,4);
        reg3 = CC3.PixelIdxList; % regions
        for i = 1:numel(reg3)
            reg4 = reg3{i};
            rc2 = zeros(size(Rc));
            rc2(reg4)=1;
            if ~isempty(find(pix==find(rc2>0)))
                break;            
            end
        end
        RCc = im2bw(conv2(rc2,ones(3),'same'),0);
        RC = rc2 + double(uint16(RCc - Lt));
        RCC = conv2(RC,ones(2),'same');
        RCC(RCC<4)=0;
        rcDown = circshift(RCC,1,1);
        rcRight = circshift(RCC,1,2);
        rcDownRight = circshift(rcRight,1,1);   
        RCC = im2bw(RCC+rcRight+rcDown+rcDownRight,0);
        Lt = double(uint16(Lt - RCC));
        [Bt,Lt0]=bwboundaries(Lt,'noholes');
        boundary = Bt{1};
        %figure(99); imagesc(Lt0); hold on; plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 2); hold off;
        if 1 % join structures
            Lt00 = Lt0;
            ns = numel(Bt);
            excRc =  zeros(size(Rc)); % excluded pixels
            for i = 1:ns
                szSt(i) = size(Bt{i},1);
                Lt0_ = Lt0;
                Lt0_(Lt0~=i) = 0; 
                excRc(:,:,i) = double(im2bw(Lt0_,0));
            end
            [~,mxSt] = max(szSt);
            Lt0_ = Lt0;
            Lt0_(Lt0~=mxSt) = 0; 
            LtMain = double(im2bw(Lt0_,0));
            LtMass = conv2(LtMain,ones(10),'same');
            Lt0 = double(im2bw(Lt0,0));
            for i = 1:ns % join structures
                if i == mxSt, continue; end;
                ps = [0 1 0; 1 0 1; 0 1 0]; % '+' conv window (plus sign)

                excRcCv3 = conv2(excRc(:,:,i),ones(3),'same');
                excRcCv5 = conv2(excRc(:,:,i),ones(5),'same');
                excRcCvOverLap = (excRcCv5 + 50*LtMain);
                excRcCvOverLap(excRcCvOverLap<51)=0;
                excRcCvOverLapPsCv = conv2(double(im2bw(excRcCvOverLap,0)),ps,'same');
                excRcConnPix = excRcCv3 + excRcCvOverLapPsCv*50; % connecting pixel
                excRcConnPix(excRcConnPix<51)=0;
                if numel(find(excRcConnPix>0))>1
                    excRcConnPixW = LtMass.*excRcConnPix; %weighted
                    excRcConnPix(excRcConnPixW~=max(excRcConnPixW(:)))=0;
                    Lt00(Lt00==i) = 0;
                end

                Lt0 = Lt0 + double(im2bw(excRcConnPix,0));
            end
        end

        figure(999); imagesc(Lt0); hold on; plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 2); hold off;
        [Bt,Lt02]=bwboundaries(Lt0,'noholes');

        if numel(Bt)>1 % remove small recs
            for i = 1:numel(Bt)
                ss(i) = size(Bt{i},1);
            end
            [~,mxss] = max(ss);
            Bt = Bt(mxss);
            Lt0_ = Lt02;
            Lt02 = zeros(size(Rc));
            Lt02(Lt0_==mxss) = 1;
        end

        Btn{ixs} = Bt{1};
        Ltn(:,:,ixs) = Lt02;
        if 0 % debug
            plot(boundary(:,2), boundary(:,1), 'k', 'LineWidth', 2)
            %figure(11); imagesc(Lt0); axis image
            %figure(12); imagesc(Lt); axis image

            %figure(16);imagesc(Lt(:,:,k))
            figure(19); plotBoundaries(Bt,sum(Lt,3).*R);
            ccc=3;
        end           
    end
end