% filtering
% IMG --> IMGfilt
IMGfilt2 = nan(size(IMG));
if inFocus_w_1X
    IMGfilt2_ = imfilter(IMG,sqr2x2,'symmetric');
    IMGfilt2 = imfilter(IMGfilt2_,lap5x5,'symmetric');
end
IMGfilt1_ = imfilter(IMG,gausHat1,'symmetric');
IMGfilt1 = imfilter(IMGfilt1_,lap,'symmetric');
IMGfilt = max([IMGfilt1(:) IMGfilt2(:)],[],2);
IMGfilt = reshape(IMGfilt,size(IMGfilt1));
