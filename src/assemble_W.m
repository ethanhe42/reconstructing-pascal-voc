function W = assemble_W(corr, corr_mirror, in_class, kp_ids)
    kp_obj = corr.obj_keypoints(in_class);
    
    n_obj = numel(kp_obj);
    n_kp = size(kp_obj{1},1);

    kp = zeros(2,n_kp,n_obj);
    for i=1:n_obj
        tmp = kp_obj{i}';
        thenan = tmp(1,:) == 0 & tmp(2,:) == 0;
        tmp(:,thenan) = NaN;
        kp(:,:,i) = tmp(1:2,:);
    end
    
    kp = kp(:,kp_ids,:);
    
    kp_obj_mirror = corr_mirror.obj_keypoints(in_class);
    kp_mirror = zeros(2,n_kp,n_obj);
    for i=1:n_obj
        tmp = kp_obj_mirror{i}';
        thenan = tmp(1,:) == 0 & tmp(2,:) == 0;
        tmp(:,thenan) = NaN;
        kp_mirror(:,:,i) = tmp(1:2,:);
    end     
     
    kp_mirror = kp_mirror(:,kp_ids,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Compute 3d model and cameras %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % marques & costeira sfm
    W = [];
    for j=1:n_obj
        duh = [kp(1,:,j); kp(2,:,j)];
        W = [W; duh];
    end    
    
    
    for j=1:n_obj
        duh = [kp_mirror(1,:,j); kp_mirror(2,:,j)];
        W = [W; duh];
    end
end