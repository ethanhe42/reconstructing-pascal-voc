function [obj_ids_in_class, filename_sfm, filename_ref, sel_kp_ids] = compute_and_cache_sfm_model(exp_dir, sfm_models_folder, sel_class, n_iter_1, n_iter_2, largest_rigid, imgset_all)    
    DefaultVal('*largest_rigid', 'false'); % use just keypoints in largest rigid part, instead of on whole object
    DefaultVal('*imgset_all', '''all_gt_segm_kp''');
    
    load('voc_kp_metadata.mat');
    
    mask_type = 'ground_truth';

    sb = SegmBrowser(exp_dir, mask_type, imgset_all);
    corr = CorrespBRKL(exp_dir, imgset_all, sb);        
    
    % load keypoints
    y = sb.get_overlaps_wholes(1:numel(sb.whole_2_img_ids));
    [q, class] = max(y, [], 2);
    obj_ids_in_class = find(class == sel_class);

    [trunc, occl] = sb.are_wholes_trunc(obj_ids_in_class);
    n_kp = size(corr.obj_keypoints{obj_ids_in_class(1)},1);
    
    if(numel(obj_ids_in_class(~(trunc|occl))) < 3*n_kp)
        NO_OCCL = false; % too few unnoccluded examples, use occluded ones in rigid SFM too
    else
        NO_OCCL = true;
        obj_ids_in_class(trunc | occl) = [];
    end
    
    kp_obj = corr.obj_keypoints(obj_ids_in_class);

    kp_names = metadata.kp_names{sel_class};

    n_obj = numel(kp_obj);
   
    kp = zeros(2,n_kp,n_obj);
    for i=1:n_obj
        tmp = kp_obj{i}';
        thenan = tmp(1,:) == 0 & tmp(2,:) == 0;
        tmp(:,thenan) = NaN;
        kp(:,:,i) = tmp(1:2,:);
    end
    
    if(strcmp(imgset_all, 'all_gt_segm_kp'))
        imgset_all_mirror = 'all_gt_segm_mirror_kp';
    else
        imgset_all_mirror = [imgset_all '_mirror'];
    end    
    
    sb_mirror = SegmBrowser(exp_dir, mask_type, imgset_all_mirror);
    corr_mirror = CorrespBRKL(exp_dir, imgset_all_mirror, sb_mirror);
    kp_obj_mirror = corr_mirror.obj_keypoints(obj_ids_in_class);
    kp_mirror = zeros(2,n_kp,n_obj);
    for i=1:n_obj
        tmp = kp_obj_mirror{i}';
        thenan = tmp(1,:) == 0 & tmp(2,:) == 0;
        tmp(:,thenan) = NaN;
        kp_mirror(:,:,i) = tmp(1:2,:);
    end        
    n_obj = 2*n_obj;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Compute 3d model and cameras %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp('Starting rigid SFM');    
    if largest_rigid
        [n_kp_srt,kp_srt_ids] = sort(cellfun(@(a) size(a,2), metadata.rigid_parts{sel_class}), 'descend');
        
        if 0
            % identify largest rigid component
            kp_ids = metadata.rigid_parts{sel_class}{kp_srt_ids(1)};
        else
            % identify set of largest rigid components such that number of
            % keypoints is larger than a minimum (8)
            kp_ids = [];
            counter = 1;
            while numel(kp_ids) < 8                
                kp_ids = [kp_ids metadata.rigid_parts{sel_class}{kp_srt_ids(counter)}];
                counter = counter + 1;
            end
            kp_ids = unique(kp_ids);
        end
    else
        kp_ids = 1:size(kp,2);
    end
            
    kp_considered = kp_names(kp_ids)
    
    % marques & costeira sfm
    W = assemble_W(corr, corr_mirror, obj_ids_in_class, kp_ids);
 
    kp = kp(:,kp_ids,:);
    kp_mirror = kp_mirror(:,kp_ids,:);
    
    if(~exist(sfm_models_folder, 'dir'))
        mkdir(sfm_models_folder);
    end  

    if(strcmp(imgset_all, 'all_gt_segm_kp'))
        filename_sfm = [sfm_models_folder 'saved_mc_class_' int2str(sel_class) '_it_' ...
            int2str(n_iter_1) '_' int2str(n_iter_2) '_sfm'];
    else
        filename_sfm = [sfm_models_folder 'saved_mc_class_' int2str(sel_class) '_it_' ...
            int2str(n_iter_1) '_' int2str(n_iter_2) '_sfm_' imgset_all];
    end        

    if largest_rigid
        filename_sfm = [filename_sfm '_largest_rigid'];
    end
       
    if(~exist([filename_sfm '.mat'], 'file'))
        t = tic();       
        
        if(sel_class == 5) % bottle's keypoints all lie in a plane
            [Motion, Shape, Transl, errors] = factorization_plane(W,n_iter_1,n_iter_2);
        else
            [Motion, Shape, Transl, errors] = factorization_custom(W,n_iter_1,n_iter_2);
        end
        toc(t)

        T = reshape(Transl, 2, n_obj);
        M = reshape(Motion', 6, n_obj);

        T = [T; zeros(1, size(T,2))];

        R = zeros(3,3,n_obj);
        for i=1:n_obj
            thisR = [M(1:3,i)'; M(4:6,i)'];

            % reconstruct 3rd line of  matrix
            scale = norm(thisR(1,:));
            R_3 = cross((1/scale)*thisR(1,:), (1/scale)*thisR(2,:));
            thisR = [thisR; scale*R_3];

            R(:,:,i) = thisR;
        end

        [R,Shape] = flip_shape(R,Shape,cat(3,kp, kp_mirror),metadata.right_coordinate_sys{sel_class},kp_ids);
        save([filename_sfm '.mat'], 'Motion', 'Shape', 'Transl', 'T', 'R', 'errors', 'kp_ids');
    else
        load([filename_sfm '.mat']);
    end    
    sel_kp_ids = kp_ids;
        
    if ~NO_OCCL
        [trunc, occl] = sb.are_wholes_trunc(obj_ids_in_class);
        obj_ids_in_class(trunc | occl) = [];
        T(:,[trunc|occl; trunc|occl]) = [];
        R(:,:,[trunc|occl; trunc|occl]) = [];    
        kp(:,:,trunc|occl) = [];
        kp_mirror(:,:,trunc|occl) = [];
    end
    
    if 0 
        show_3d_model(Shape, kp_names, 'convex_hull'); axis equal;
    end

    if 0
        obj_id = 4;
        sb.show_wholes(obj_ids_in_class(obj_id)); hold on;
        duh = bsxfun(@plus, R(:,:,obj_id)*Shape, T(:,obj_id));
        cmap = jet(size(duh,2));
        for k=1:size(duh,2)
            hold on; plot3(duh(1,k), duh(2,k), duh(3,k),  'o', 'MarkerSize', 18, 'LineWidth', 10, 'Color', cmap(k,:));
        end
    end
    
    disp('Finished rigid SFM, refining...');
    REFINE_CAMERAS = true;
    if REFINE_CAMERAS
        %
        % refine cameras
        %
        filename_ref = [filename_sfm '_refine'];
        if(~exist([filename_ref '.mat'], 'file'))
            all_masks = sb.get_masks(obj_ids_in_class);

            t = tic();
            newR = zeros(size(R,1), size(R,2), size(kp,3));
            newT = zeros(size(T,1), size(kp,3));
            parfor i=1:size(kp,3)
                try
                    [newR(:,:,i),newT(:,i)] = refine_all_cameras(Shape,R(:,:,i),T(:,i), all_masks(i), kp(:,:,i));
                catch
                    disp('failed on one camera');
                end                
            end
            toc(t)

            all_masks_mirror = sb_mirror.get_masks(obj_ids_in_class);

            t = tic();      
            newR_mirror = zeros(size(R,1), size(R,2), size(kp,3));
            newT_mirror = zeros(size(T,1), size(kp,3));
            c = size(kp,3);
            parfor i=(c+1):size(R,3)
                try
                    [newR_mirror(:,:,i-c),newT_mirror(:,i-c)] = refine_all_cameras(Shape,R(:,:,i),T(:,i), all_masks_mirror(i-c), kp_mirror(:,:,i-c));
                catch
                    disp('failed on one camera');
                end
            end
            toc(t)

            R = cat(3, newR, newR_mirror);
            T = cat(2, newT, newT_mirror);

            save([filename_ref '.mat'], 'R', 'T');   

            disp('Camera refinement finished.');
        end
    end
end
    
