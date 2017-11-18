    % dimension of of the oversampled super-resolution image
    debug_mode =0;


    
    if(debug_mode)
        mag_full = 2 / x_dim ;
        min_contrib = 1; % min CS recovered value to be considered as a successful molecule identification
    end
    
    x_dim_est = x_dim * div;
    y_dim_est = y_dim * div;

    % generate global grid for final image estimate
    % these dimensions should correspond to your raw and final image
    x_pos_est_full = linspace(-1,1,x_dim_est);
    y_pos_est_full = linspace(-1,1,y_dim_est);
    x_inx_full = x_pos_est_full(x_dim_est/x_dim/2 : x_dim_est/x_dim : end);
    y_inx_full = y_pos_est_full(y_dim_est/y_dim/2 : y_dim_est/y_dim : end);

    [x_inx_full y_inx_full] = meshgrid(x_inx_full, y_inx_full);
    [x_pos_est_full y_pos_est_full] = meshgrid(x_pos_est_full, y_pos_est_full);

    % generate grid for local estimation
    x_pos_est = linspace(-1,1, div*boxsize);
    y_pos_est = linspace(-1,1, div*boxsize);

    x_inx = x_pos_est((div/2) : div : end);
    y_inx = y_pos_est((div/2) : div : end);

    % provide pixel padding for local CS estimate
    % This should correspond to the box size around peak intensity
    % -- Begin padding the estimated local patch
    dx__ = x_pos_est(2)-x_pos_est(1);
    dy__ = y_pos_est(2)-y_pos_est(1);

    for mm = 1:margin
        x_pos_est = [x_pos_est(1)-dx__ x_pos_est];
        x_pos_est = [x_pos_est x_pos_est(end)+dx__];
        y_pos_est = [y_pos_est(1)-dy__ y_pos_est];
        y_pos_est = [y_pos_est y_pos_est(end)+dy__];
    end

    [x_inx y_inx] = meshgrid(x_inx,y_inx);
    [x_pos_est y_pos_est] = meshgrid(x_pos_est,y_pos_est);

    % generate measurement matrix
    len = length(x_pos_est(:));
    %A = zeros(x_inx, len + 1);
    for ii = 1:len
        img_kernel = MolKernel(x_inx, y_inx, x_pos_est(ii), y_pos_est(ii), mag);
        A(:,ii) = img_kernel(:);
    end
    
    c__ = sum(A);
    PSF_integ = max(c__); % integration of the PSF over space, used for normalization
    c__ = c__./ PSF_integ; % normalize to 1
    A = A./ PSF_integ;
    
    % add the extra optimization variable for the estimation of the background
    len = len + 1;
    c__(len) = 0;
    A(:,len) = 1;
    
    len = size(A,2);
    A = sparse(A);
    
    % begin to formulate optimization
    cvx_quiet(true);
    
    