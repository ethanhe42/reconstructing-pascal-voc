function reconstruct_pascal_getall(sel_class, n_iter_1, n_iter_2, N_SAMPLES_PER_OBJ, MAX_DEV, name, refinement, imprinting)
    DefaultVal('*refinement', 'true');           
    DefaultVal('*imprinting', 'true');    
    
    exp_dir = add_all_paths();

    mask_type = 'ground_truth';
    imgset_pascal = 'all_gt_segm_kp';
    imgset_pascal_mirror = 'all_gt_segm_mirror_kp'; 

    imgset_all = imgset_pascal;
    imgset_all_mirror = imgset_pascal_mirror;
    
    load('voc_kp_metadata.mat');

    kp_names = metadata.kp_names{sel_class};
   
    sb_all = SegmBrowser(exp_dir, mask_type, imgset_all);
    sb_all_mirror = SegmBrowser(exp_dir, mask_type, imgset_all_mirror);

    sfm_folder = ['./Results/'];

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
        show_cameras_sv(Shape, kp_names(sel_kp_ids), R(:,:,viz_ids), T(:,viz_ids), Is, the_masks);
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

    % collect masks, should perhaps fill any holes
    all_masks = sb_all.get_masks(in_class);
    all_sym_masks = sb_all_mirror.get_masks(in_class);                
    
    really_all_masks = [all_masks; all_sym_masks];
    
    if 0
        areas = compute_mask_areas(sb_all);
        areas = areas(in_class);
    end
           
    max_angle_deviations = [MAX_DEV MAX_DEV MAX_DEV];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Get data of images having ground truth meshes %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    test_mask = 'ground_truth';
    sb_test = SegmBrowser(exp_dir, test_mask, imgset_pascal);
    corr_test = CorrespBRKL(exp_dir, imgset_pascal, sb_test);

    y = sb_test.get_overlaps_wholes(1:numel(sb_test.whole_2_img_ids));
    [q, class] = max(y, [], 2);
    in_class_test = find(class == sel_class);
        
    sb_test_mirror = SegmBrowser(exp_dir, test_mask, imgset_pascal_mirror);
    corr_test_mirror = CorrespBRKL(exp_dir, imgset_pascal_mirror, sb_test_mirror);
    
    kp_test = corr_test.obj_keypoints(in_class_test);
    kp_test = cellfun(@(a) a(:,1:2), kp_test, 'UniformOutput', false);    
    kp_test_mirror = corr_test_mirror.obj_keypoints(in_class_test);
    kp_test_mirror = cellfun(@(a) a(:,1:2), kp_test_mirror, 'UniformOutput', false);
    
    all_kp_test = [kp_test kp_test_mirror];
    masks_test = sb_test.get_masks(in_class_test);
    masks_test_mirror = sb_test_mirror.get_masks(in_class_test);
    all_masks_test = [masks_test; masks_test_mirror];
    
    % gets ids of gt within the pascal+gt set
    ids_all = sb_all.whole_2_img_ids(in_class);
    test_ids = [ids_all; ids_all+(size(R,3)/2)];    
    pascal_ids = test_ids;

    R_test = R;    
    T_test = T;
    R_pascal = R;
    T_pascal = T;
    
    kp_pascal = kp_obj;
    kp_pascal_mirror = kp_obj_mirror;
    
    all_kp_pascal = [kp_pascal kp_pascal_mirror];
    masks_pascal = sb_all.get_masks(in_class);
    masks_pascal_mirror = sb_all_mirror.get_masks(in_class);
    all_masks_pascal = [masks_pascal; masks_pascal_mirror];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Sample triplets of object images %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Reconstruct object using a single view from the test
    % object and 2 others from a pool of examples of the same category from
    % the pascal voc data
    [test_ids_sel, pascal_ids_sel, axs] = view_comb_sampling_general(R_test, R_pascal, N_SAMPLES_PER_OBJ, max_angle_deviations);     
    axis_masks = get_average_mask(axs, R, T, really_all_masks);
        
    ids = [test_ids_sel pascal_ids_sel];
    
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
                
    out_folder = ['SFM_models/' name '/' sb_all.categories{sel_class} '/'];
    mkdir(out_folder);

    save([out_folder 'reconstruction_data'], 'in_class', 'imgset_all', 'ids', 'N_VOXELS', 'max_angle_deviations','-v7.3');
    obj_ids = unique(test_ids_sel);
    for i=1:numel(obj_ids)
        tic
        disp(i);
        file_name = [out_folder int2str(i) '.mat'];
        if(exist(file_name,'file'));
           continue
        end
        
        ids_sel = ids(ids(:,1)==obj_ids(i),:);                
        ids_sel = unique(ids_sel, 'rows');
                                                         
        reconstructions = sel_best_reconstruction_getall(ids_sel,R,T,all_kp,really_all_masks, N_VOXELS,axis_masks,flags);
        
        if 0
            t = reconstructions.triples(1:end/2);
            
            t(t==0) = [];
            
            the_masks = [all_masks(t); all_sym_masks(t)];
            
            figure;
            Is = [sb_all.get_Imgs(sb_all.whole_2_img_ids(in_class(t)));sb_all_mirror.get_Imgs(sb_all_mirror.whole_2_img_ids(in_class(t)))];
            subplot_auto_transparent_multiple_imgs(the_masks, Is); 
        end

        if(numel(reconstructions) == 1 && isempty(reconstructions.triples))           
            theR = [];
            theT = [];
        else
            allids = arrayfun(@(a) a.triples(1), reconstructions);
            assert(numel(unique(allids))==1);
            theR = R(:,:,allids(1));
            theT = T(:,allids(1));

            for j=1:numel(reconstructions)
                reconstructions(j).vertices = bsxfun(@plus,theR*reconstructions(j).vertices',theT);
            end
        end
        
        if 0
            I = sb_all.get_Imgs(sb_all.whole_2_img_ids(in_class(ids_sel(1,1))));
            I = I{1};
            
            % visualize reconstruction 
            sc(I); hold on;
            
            rec_id = 10;
            trisurf(reconstructions(rec_id).faces, reconstructions(rec_id).vertices(1,:), ...
                reconstructions(rec_id).vertices(2,:), reconstructions(rec_id).vertices(3,:), 'FaceColor', 'blue');
            axis equal;
            set(gca,'zdir','reverse');
        end
        save(file_name, 'reconstructions','ids_sel', 'theR', 'theT', '-v7.3');
        toc

        close all;
    end
    time = toc(t_total)
    save([out_folder 'reconstruction_data'], 'in_class', 'imgset_all', 'N_VOXELS', 'max_angle_deviations','time','flags', 'Shape', 'axis_masks', '-v7.3');
end

 
function [Xvol, Yvol, Zvol, dt] = precompute_data(R, T, Shape, really_all_masks, step, border, faster)
    %TODO: this should probably be adjusted depending on the class, to ensure
    %that the volumetric representation covers all the volume it needs.    
    mi = floor(min(Shape,[],2)) -border;
    ma = ceil(max(Shape,[],2)) + border;        
    [Xvol,Yvol,Zvol] = meshgrid(mi(1):step:ma(1),mi(2):step:ma(2),mi(3):step:ma(3));    
    vol = [Xvol(:)'; Yvol(:)'; Zvol(:)'];
    dt = zeros(size(vol,2),size(R,3));
    
    %parfor h=1:size(R,3)
    for h=1:size(R,3)
        thisR = R(:,:,h);
        thisT = T(:,h);
        mask = really_all_masks{h};
        
        rotVol = thisR*vol;
        rotVol =  bsxfun(@plus,rotVol,thisT);
        dt(:,h) = distance_from_mask(mask,rotVol,faster);
    end
end

function dt = distance_from_mask(mask,rotVol,faster)
% Returns a negative value for voxels in the volume and positive value
% outside.


[M,N] = size(mask);

minX = floor(min(rotVol(1,:))-2);
minY = floor(min(rotVol(2,:))-2);

maxX = ceil(max(rotVol(1,:))+2);
maxY = ceil(max(rotVol(2,:))+2);

cX = max(0,-minX+1);
cY = max(0,-minY+1);

aux_mask = false(cY + max(maxY,M),cX + max(maxX,N));

aux_mask(cY+1:cY+M, cX+1:cX+N) = mask;

rotVol(1,:) = rotVol(1,:) + cX;
rotVol(2,:) = rotVol(2,:) + cY;

ind = sub2ind(size(aux_mask),round(rotVol(2,:)),round(rotVol(1,:)));

if(faster)
    dt = -2*double(aux_mask(ind))' + 1;
else
    dt = zeros(size(rotVol,2),1);

    [x,y] = getXYfromMask(aux_mask);
    
    %Discrete approximation of distance, it simply takes the closest point in
    %terms of the 2D distance transform, and then refines that distance.
    [~,IDX] = bwdist(aux_mask);
    sel_points = IDX(ind);
    d = sqrt((x(sel_points)-rotVol(1,:)').^2 + (y(sel_points)-rotVol(2,:)').^2);
    
    dt(~aux_mask(ind)) = d(~aux_mask(ind));
    
    [~,IDX] = bwdist(~aux_mask);
    sel_points = IDX(ind);
    d = sqrt((x(sel_points)-rotVol(1,:)').^2 + (y(sel_points)-rotVol(2,:)').^2);
    
    dt(aux_mask(ind)) = -d(aux_mask(ind));
end

end


function [x,y,m] = getXYfromMask(mask)

[x,y] = meshgrid(1:size(mask,2),1:size(mask,1));

x = x(:);
y = y(:);
m = mask(:);


end

function out = assig2cell(assig)
    un = unique(assig);
    n = numel(un);
    
    out = cell(n,1);
    for i=1:n
        out{i} = find(assig==un(i));
    end
    out(isempty(out)) = [];
end
