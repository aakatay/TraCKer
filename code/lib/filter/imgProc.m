close 
clear
load('J')
load('Coeff')


        gaus=fspecial('gaussian', 5, 1);
        lap=[-1,-1,-1;-1,8,-1;-1,-1,-1];
        A = 1;  
        parab = A*[sqrt(2), 1, sqrt(2); 1, 0, 1; sqrt(2), 1, sqrt(2) ];
        parabFilt = A*sqrt(2) - parab;
%% raw image        
figure
subplot(3,5,1)
imagesc(J); axis image;
title('raw data')
[Boy1,En1]=size(J);

%% low pass
Norder = 6; cutFreq = 0.5;
b = fir1(Norder,cutFreq);h = b' * b;
J2 = imfilter(J, h);
subplot(3,5,2)
imagesc(J2); axis image;
title('low pass')

med = median(J2(:));
J3=1.0*(J2-med/1)/1;
J3=J2;

%% back ground
%J3 = imfilter(J2,lap,'symmetric');
subplot(3,5,3)
imagesc(J3); 
axis image;
title('laplacian')

%% regional max & boundary detection

C1 = 550; C2 = 900;
dCoeff = (C2-C1)/4;
coeffVec = C1:dCoeff:C2;

pIx = 1;
for Coeff = coeffVec
    N = J3/Coeff; 
    imgMax = max(N(:));
    N = N.^2/imgMax; clear BW;
    BW = imregionalmax(N, 8);
    subplot(3,5,5+pIx)
    imagesc(BW); axis image;
    title('regional max')

    %BW = J3;
    clear B L;
    [B,L] = bwboundaries(BW,'noholes');

    for m=1:length(B) % for each spot
        c=cell2mat(B(m));
        %M(size(J),c(:,1),c(:,2)) = 1;
        %csize=(max(c(:,1))-min(c(:,1)))*(max(c(:,2))-min(c(:,2)));
        Py(m)=uint16(mean(c(:,1)));
        Px(m)=uint16(mean(c(:,2)));
    end
    
    % pits
    subplot(3,5,5*2+pIx)
    scatter(Px,Boy1-Py,'.');axis image;
    xlim([1 En1]);
    ylim([1 Boy1]);
    title('pits')
    
    pIx = pIx + 1;
end

%% pits

subplot(3,5,5)
scatter(Px,Boy1-Py,'.');axis image;
xlim([1 En1]);
ylim([1 Boy1]);
title('pits')
maximize