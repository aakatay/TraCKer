
lowresimg = (IMG - ccd_base_line) * photon_per_count;

if debug_mode
   lowresimg_display = lowresimg;
end
img_recover = zeros(x_dim_est, y_dim_est);

%scan the optimization box throughout image
pre_xloc = 1 + scanoverlap - boxsize;
pre_yloc = 1;
while true;
   pre_xloc = pre_xloc + boxsize - scanoverlap;
   if pre_xloc + boxsize - 1 > x_dim
       pre_yloc = pre_yloc + boxsize - scanoverlap;
       pre_xloc = 1;
       if pre_yloc + boxsize - 1 > y_dim 
           break;
       end
   end

   % The boundary of the optimization box
   inx_x_l = pre_xloc;
   inx_x_u = pre_xloc + boxsize - 1;
   inx_y_l = pre_yloc;
   inx_y_u = pre_yloc + boxsize - 1;

   if debug_mode
      imagesc((lowresimg_display)), colormap(gray), hold on;
      rectangle('Position',[inx_x_l-0.5,inx_y_l-0.5,boxsize,boxsize],'linewidth',1,'EdgeColor','y'),pause(0.5);
   end

   % crop from the original image
   cropraw = lowresimg(inx_y_l:inx_y_u, inx_x_l:inx_x_u);

   L1sum = norm(cropraw(:),1);
   L2sum = sqrt(L1sum); % Target of least square estimation based on Poisson statistics

   %disp(sprintf('a1:%.02f a2:%.02f',a1,a2))
   
    % checks if std sum is smaller than total intensity squared
   % multipl by optim. coeff
   if norm(cropraw(:)-L1sum/length(cropraw(:)),2) < L2sum * red_chi2;
       continue;    % no molecules in the patch, just flat background
   end

   % optimization using CVX
   %Weighting = diag(1./sqrt(abs(cropraw(:))+1));
   %img_est = MolEst_eps(Weighting*cropraw(:),len,Weighting*A,c__,(boxsize*boxsize-1)*red_chi2);
   img_est = MolEst_eps(cropraw(:),len,A,c__,L2sum * red_chi2);
   est_background = img_est(len);
   if(debug_mode)
       fprintf('BG=%.1f, L1=%.1f (target: %.1f), L2=%.1f (target: %.1f), ',...
              est_background,...
              norm(c__*img_est,1),             L1sum-(boxsize*boxsize)*est_background,...
              norm(A*img_est-cropraw(:),2),  L2sum);
   end

   img_est = reshape(img_est(1:len-1), [length(x_pos_est) length(y_pos_est)]);

   % remove margin
   img_est = img_est(margin+1:end-margin, margin+1:end-margin);
   if extramargin > 0
      img_est([1:extramargin*div end-extramargin*div+1:end],:) = 0;
      img_est(:,[1:extramargin*div end-extramargin*div+1:end]) = 0;
   end

   % put this hi-res patch of image back into collection array
   hr_x_l = (inx_x_l-1)*div+1;
   hr_y_l = (inx_y_l-1)*div+1;
   hr_x_u = hr_x_l + (div*boxsize)-1;
   hr_y_u = hr_y_l + (div*boxsize)-1;

   img_recover_temp = zeros(x_dim_est,y_dim_est);
   img_recover_temp(hr_y_l:hr_y_u, hr_x_l:hr_x_u) = img_est;

   % add the optimized patch to the recovered super-resolution image
   img_recover = img_recover + img_recover_temp;

   % display othe optimization procedure for debugging only
   if debug_mode
%               inx_contribute = find(img_recover_temp(:) > min_contrib);
       inx_contribute = find(img_recover_temp(:) >= min_contrib);
       if isempty(inx_contribute)
           fprintf('nothing found.\n');
           continue;
       else
           fprintf('%d contribution points.\n',length(inx_contribute(:)));
       end

       % calculate the recovered low resolution image based on the optimization results
       lowresimg_temp = zeros(size(IMG));
       for jj = 1:length(inx_contribute)
          h__ = MolKernel(x_inx_full, y_inx_full, ...
                    x_pos_est_full(inx_contribute(jj)),...
                    y_pos_est_full(inx_contribute(jj)), mag_full);
          lowresimg_temp = lowresimg_temp + img_recover_temp(inx_contribute(jj))*h__;
       end

       lowresimg_display = lowresimg_display - lowresimg_temp / PSF_integ;

       imagesc(lowresimg_display); drawnow;
   end
end


% out : img_recover
        