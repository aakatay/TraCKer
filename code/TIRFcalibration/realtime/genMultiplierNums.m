clear all;
close all;

    sn  = 3; % set number
    FN  = sprintf('set%i.mat',sn);
    
    mNlim = [2000 4000 6000 8000];
    np = 2000;
    
    mns(1,:) = mNlim(1)+1:mNlim(1)+np;
    % set 2
    mns(2,:) = mNlim(2)+1:mNlim(2)+np;
    % set 3
    mns(3,:) = mNlim(3)+1:mNlim(3)+np;

    mnsMul2D(:,:,1) = mns(1,:)' * mns(2,:); % pns12
    pns12 = mnsMul2D(:,:,1);
    mnsMul(1,:) = pns12(:)'; 
    
    mnsMul2D(:,:,2) = mns(2,:)' * mns(3,:); % pns23
    pns23 = mnsMul2D(:,:,2);
    mnsMul(2,:) = pns23(:)'; 
    
    mnsMul2D(:,:,3) = mns(3,:)' * mns(1,:); % pns31
    pns31 = mnsMul2D(:,:,3);
    mnsMul(3,:) = pns31(:)'; 
    
    sIx = [1 2 ; 2 3 ; 3 1]; % set index
    
    iFrst = 1;
    if exist(FN)
        load(FN); 
        ilast = numel(fct1);
        iFrst = ilast+1;
    end
    
    hw = waitbar(0);
    for i = iFrst:np^2
        [y,x]=find(mnsMul2D(:,:,sn)==mnsMul(sn,i));
        if numel(y)>1
            diff = y-x;
            dix = diff==0;
            fdix = find(dix);
            if ~isempty(fdix)
                ix = fdix;
                fct1(i) = mns(sIx(sn,1),y(ix));
                fct2(i) = mns(sIx(sn,2),x(ix));
            else
                yx = y.*x;
                [~,ix] = max(yx); % select indices close to diagonal
                fct1(i) = mns(sIx(sn,1),y(ix));
                fct2(i) = mns(sIx(sn,2),x(ix));
            end
            
            ccc=3;
        else
            fct1(i) = mns(sIx(sn,1),y);
            fct2(i) = mns(sIx(sn,2),x);
        end
        if rem(i,1000)==0
            wtxt = sprintf('%i/%i',i/1000,np^2/1000);
            waitbar(i*1/np^2,hw,wtxt)
        end
    end
    close(hw)
    save(FN,'fct1','fct2')