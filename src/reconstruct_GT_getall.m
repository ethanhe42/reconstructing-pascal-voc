function reconstruct_GT_getall(sel_class, n_iter_1, n_iter_2, N_SAMPLES_PER_OBJ, MAX_DEV, name, refinement, imprinting)        
    exp_dir = add_all_paths();

    DefaultVal('*refinement', 'true');           
    DefaultVal('*imprinting', 'true');
    
    mask_type = 'ground_truth';

    imgset_all = 'all_gt_segm_kp_synth_5_views_per_mesh';
    imgset_all_mirror = 'all_gt_segm_kp_synth_5_views_per_mesh_mirror';

    load('voc_kp_metadata.mat');

    kp_names = metadata.kp_names{sel_class};
   
    sb_all = SegmBrowser(exp_dir, mask_type, imgset_all);
    sb_all_mirror = SegmBrowser(exp_dir, mask_type, imgset_all_mirror);

    %%%%% This select class images, initializes cameras for all objects using rigid SFM and refines them using the silhouettes %%%            
    sfm_folder = ['./SFM_models/'];

    if(metadata.articulated(sel_class))
        largest_rigid = true;        
    else
        largest_rigid = false;
    end

    [in_class, filename_sfm, filename_ref, sel_kp_ids] = compute_and_cache_sfm_model(exp_dir, sfm_folder, sel_class, n_iter_1, n_iter_2, largest_rigid, imgset_all);            
    load([filename_sfm '.mat']);
    if(refinement)
        var = load([filename_ref '.mat'], 'R', 'T');        
        R = var.R;
        T = var.T;    
    end
    
    [Shape, R, k] = postprocess_sfm(Shape, R);

    if 0
        show_3d_model(Shape, kp_names(sel_kp_ids), 'convex_hull'); axis equal;
    end

    if 0
        % visualize some cameras to see if it's ok
        figure;
        viz_ids = 1:20;        
        Is = sb_all.get_Imgs(sb_all.whole_2_img_ids(in_class(viz_ids)));
        the_masks = sb_all.get_masks(in_class(viz_ids));
        show_cameras(Shape, kp_names(sel_kp_ids), R(:,:,viz_ids), T(:,viz_ids), Is, the_masks);
    end

    corr = CorrespBRKL(exp_dir, imgset_all, sb_all);
    corr_mirror = CorrespBRKL(exp_dir, imgset_all_mirror, sb_all_mirror);
    kp_obj = corr.obj_keypoints(in_class);
    kp_obj = cellfun(@(a) a(:,1:2), kp_obj, 'UniformOutput', false);    

    kp_obj_mirror = corr_mirror.obj_keypoints(in_class);
    kp_obj_mirror = cellfun(@(a) a(:,1:2), kp_obj_mirror, 'UniformOutput', false);
    
    all_kp = [kp_obj kp_obj_mirror];
    for i=1:numel(all_kp)
        all_kp{i}(isnan(all_kp{i})) = 0;
    end
    
    n_kp = cellfun(@(a) sum(a(:,1)~=0 | a(:,2) ~= 0), kp_obj)';

    all_masks = sb_all.get_masks(in_class);
    all_sym_masks = sb_all_mirror.get_masks(in_class);                
    
    really_all_masks = [all_masks; all_sym_masks];
    
    if 0
        areas = compute_mask_areas(sb_all);
        areas = areas(in_class);
    end
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Get data of images having ground truth meshes %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    synth_mask = 'ground_truth';
    sb_synth = SegmBrowser(exp_dir, synth_mask, 'synth_5_views_per_mesh');
    corr_synth = CorrespBRKL(exp_dir, 'synth_5_views_per_mesh', sb_synth);

    y = sb_synth.get_overlaps_wholes(1:numel(sb_synth.whole_2_img_ids));
    [q, class] = max(y, [], 2);
    in_class_synth = find(class == sel_class);
        
    sb_synth_mirror = SegmBrowser(exp_dir, synth_mask, 'synth_5_views_per_mesh_mirror');
    corr_synth_mirror = CorrespBRKL(exp_dir, 'synth_5_views_per_mesh_mirror', sb_synth_mirror);
    
    kp_synth = corr_synth.obj_keypoints(in_class_synth);
    kp_synth = cellfun(@(a) a(:,1:2), kp_synth, 'UniformOutput', false);    
    kp_synth_mirror = corr_synth_mirror.obj_keypoints(in_class_synth);
    kp_synth_mirror = cellfun(@(a) a(:,1:2), kp_synth_mirror, 'UniformOutput', false);
    
    all_kp_synth = [kp_synth kp_synth_mirror];
    masks_synth = sb_synth.get_masks(in_class_synth);
    masks_synth_mirror = sb_synth_mirror.get_masks(in_class_synth);
    all_masks_synth = [masks_synth; masks_synth_mirror];
    
    % gets ids of gt within the pascal+gt set
    ids_all = sb_all.whole_2_img_ids(in_class);
    [~,the_ids,~] = intersect(sb_all.img_names(ids_all), sb_synth.img_names(in_class_synth));
    synth_ids = [the_ids; the_ids+(size(R,3)/2)];
    pascal_ids = [setdiff(1:size(R,3), synth_ids)];
    R_synth = R(:,:,synth_ids);    
    T_synth = T(:,synth_ids);
    R_pascal = R(:,:,pascal_ids);
    T_pascal = T(:,pascal_ids);
    
    pascal_ids_nonmirr = pascal_ids(1:(numel(pascal_ids)/2));
    kp_pascal = kp_obj(pascal_ids_nonmirr);
    kp_pascal_mirror = kp_obj_mirror(pascal_ids_nonmirr);
    
    all_kp_pascal = [kp_pascal kp_pascal_mirror];
    masks_pascal = sb_all.get_masks(in_class(pascal_ids_nonmirr));
    masks_pascal_mirror = sb_all_mirror.get_masks(in_class(pascal_ids_nonmirr));
    all_masks_pascal = [masks_pascal; masks_pascal_mirror];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Sample triplets of object images %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Reconstruct synthetic object using a single view from the synthetic
    % object and 2 others from a pool of examples of the same category from
    % the pascal voc data
    max_angle_deviations = [MAX_DEV MAX_DEV MAX_DEV];
    [synth_ids_sel, pascal_ids_sel, axs] = view_comb_sampling_general(R_synth, R_pascal, N_SAMPLES_PER_OBJ, max_angle_deviations);
    axis_masks = get_average_mask(axs, R, T, really_all_masks);
    
    % these ids are with respect to synth and non-synth
    % get triplets of global ids
    ids = [synth_ids(synth_ids_sel) pascal_ids(pascal_ids_sel)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Reconstruct object instances %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N_VOXELS = 60;

    if largest_rigid
        filename_precomp = ['./SFM_models/precomp_dt_' sb_all.categories{sel_class} '_' imgset_all '_largest_rigid.mat'];
    else
        filename_precomp = ['./SFM_models/precomp_dt_' sb_all.categories{sel_class} '_' imgset_all '.mat'];
    end
       
    t_total = tic();
           
    flags.is_articulated = metadata.articulated(sel_class);
    flags.rot_symmetry = metadata.view_based(sel_class);
    flags.imprinted = imprinting;
    
    if(flags.rot_symmetry)
        flags.rot_axis = get_rot_axis(Shape,metadata.bottom_top{sel_class});
    end
    flags.angle_step = pi/2;
                
    out_folder = ['Results/' name '/' sb_all.categories{sel_class} '/'];
    mkdir(out_folder);

    save([out_folder 'reconstruction_data'], 'in_class', 'imgset_all', 'ids', 'N_VOXELS', 'max_angle_deviations','-v7.3');
    for i=1:numel(unique(synth_ids_sel))
        tic
        disp(['Reconstructing image ' num2str(i)]);
        file_name = [out_folder int2str(i) '.mat'];
        
        ids_sel = ids(ids(:,1)==synth_ids(i),:);                
        [reconstructions,statistics] = sel_best_reconstruction_getall(ids_sel,R,T,all_kp,really_all_masks, N_VOXELS,axis_masks,flags);

        allids = arrayfun(@(a) a.triples(1), reconstructions);
        assert(numel(unique(allids))==1);
        theR = R(:,:,allids(1));
        theT = T(:,allids(1));
        
        for j=1:numel(reconstructions)
            reconstructions(j).vertices = bsxfun(@plus,theR*reconstructions(j).vertices',theT);
        end
        save(file_name, 'reconstructions','statistics','ids_sel', 'theR', 'theT', '-v7.3');
        toc

        close all;
    end
    time = toc(t_total)
    save([out_folder 'reconstruction_data'], 'in_class', 'imgset_all', 'ids', 'N_VOXELS', 'max_angle_deviations','Shape','time','axis_masks','flags','-v7.3');
end
