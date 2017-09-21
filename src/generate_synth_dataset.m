function generate_synth_dataset(exp_dir)
    load('voc_kp_metadata.mat');
    eval('config');
    
    % make it repeatable
    rng(1234);

    %%%%%%%%%%%%%%%%%%% params & directories & filenames %%%%%%%%%%%%%%%%%%        
    sfm_folder = ['./SFM_models/'];
    
    mask_type = 'ground_truth';
    imgset_name = ['synth_' int2str(N_SAMPLES_PER_MESH) '_views_per_mesh'];
    
    mask_dir = [exp_dir 'MySegmentsMat/' mask_type '/'];
    if(~exist(mask_dir, 'dir'))
        mkdir(mask_dir);
    end
    
    cam_mesh_dir = [exp_dir 'MyMeshes/' mask_type '/'];
    if(~exist(cam_mesh_dir, 'dir'))
        mkdir(cam_mesh_dir);
    end     
    
    file_fd = fopen([exp_dir 'ImageSets/Segmentation/' imgset_name '.txt'], 'w');    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:20
        categs{i} = VOC09_id_to_classname(i);
    end

    sb = SegmBrowser(exp_dir, 'ground_truth', 'all_gt_segm_kp');
    corr = CorrespBRKL(exp_dir, 'all_gt_segm_kp', sb);        
    
    data_dir = './Dataset/SynthMeshes/';

    y = sb.get_overlaps_wholes(1:numel(sb.whole_2_img_ids));
    [q, class] = max(y, [], 2);

    for h=1:numel(categs)        
        mesh_files = dir([data_dir categs{h} '/*_mesh.mat']);
        n_meshes = numel(mesh_files);
        assert(n_meshes>0);
        
        t_mesh = tic();
        for i=1:n_meshes        
            base_output_img_name = strtok(mesh_files(i).name,'_');
            for j=1:N_SAMPLES_PER_MESH
                output_img_name = [base_output_img_name '_view_' int2str(j)];
                fprintf(file_fd, '%s\n', output_img_name);
            end
        end
    end
    fclose(file_fd); % close imgset text file
  
    for h=1:numel(categs)
        kp_names = metadata.kp_names{h};
        n_kp = numel(kp_names);

        %%% retrieve cameras on pascal images        
        if(metadata.articulated(h))
            largest_rigid = true;
        else
            largest_rigid = false;
        end
        [in_class, filename_sfm, filename_ref, sel_kp_ids] = compute_and_cache_sfm_model(exp_dir, sfm_folder, h, n_iter_1, n_iter_2, largest_rigid);
        load([filename_sfm '.mat'], 'Shape');
        load([filename_ref '.mat'], 'R', 'T'); % R and T include versions for left-right mirror images, stacked at the end

        % there are some bad cameras (because of lack of annoatated keypoints, etc.) which we will ignore
        bad_cams = find(sum(abs(T),1) == 0);

        [Shape, R, k] = postprocess_sfm(Shape, R);
                
        if 0
            corr.show_wholes_and_keypoints(in_class(1:10), false);
        end
        
        Imgs = sb.get_Imgs(sb.whole_2_img_ids(in_class));
        dims = cellfun(@size, Imgs, 'UniformOutput', false);
        dims = cell2mat(dims);
        dims = dims(:,1:2);        

        mesh_files = dir([data_dir categs{h} '/*_mesh.mat']);
        n_meshes = numel(mesh_files);
        assert(n_meshes>0);

        t_mesh = tic();
        for i=1:n_meshes
            base_output_img_name = strtok(mesh_files(i).name,'_');

            if i < 10
                number_str = ['0' int2str(i)];
            else
                number_str = int2str(i);
            end
            
            mesh_filename = [data_dir categs{h} '/' categs{h} number_str '_mesh.mat'];
            load(mesh_filename, 'xyz', 'nodes');                        

            tri_rep = TriRep(nodes',xyz');

            % get the annotated 3D keypoints on the meshes
            kp_filename = [data_dir categs{h} '/' categs{h} number_str '_kp.mat'];
            var = load(kp_filename, 'kp');
            kp3d = var.kp;
            
            [avail_cams] = setdiff(1:(size(R,3)/2), bad_cams);
            
            for j=1:N_SAMPLES_PER_MESH 
                output_img_name = [base_output_img_name '_view_' int2str(j)];
                
                all_rp_cam = avail_cams(randperm(numel(avail_cams)));

                while 1 % while selected view does not fit in the image, keep selecting the next one
                    rp_cam = all_rp_cam(1);
                    all_rp_cam(1) = [];
                    
                    kp_obj = corr.obj_keypoints(in_class(rp_cam));
                    kp_obj = kp_obj{1};
                    
                    camR = R(:,:,rp_cam);
                    camT = T(:,rp_cam);
                    sel_dims = dims(rp_cam,:);
                    
                    % estimate rigid transformation between the mesh and the SFM model
                    anan = isnan(kp3d(1,sel_kp_ids));
                    
                    [D,Z,transf] = procrustes(Shape(:,~anan)', kp3d(:,sel_kp_ids(~anan))', 'Reflection', false);
                    global_R = transf.T';
                    global_s = transf.b;
                    global_T = transf.c(1,:)';                                        

                    transf_xyz = bsxfun(@plus, global_s*global_R*xyz, global_T);
                    cam_xyz = bsxfun(@plus,camR*transf_xyz,camT);

                    if(max(cam_xyz(1,:))>sel_dims(2) || min(cam_xyz(1,:))<1 || ...
                            max(cam_xyz(2,:))>sel_dims(1) || min(cam_xyz(2,:))<1)
                        % if it's not fitting inside the image quit on this
                        % camera
                        continue;
                    end

                    I = 255*ones(sel_dims(1), sel_dims(2), 3, 'uint8'); hold on;
                    f = figure;
                    set(f,'Renderer','zbuffer');
                    sc(I, 'w', true(sel_dims(1), sel_dims(2))); hold on;
                    trisurf(nodes',cam_xyz(1,:),cam_xyz(2,:),cam_xyz(3,:),'FaceColor','red');
                    set(gca,'zdir','reverse');
                    axis equal;
                    render_im = getframe(f);
                    render_im = imresize(render_im.cdata, [size(I,1) size(I,2)]);
                    assert(all(size(render_im) == size(I)));
                    imwrite(render_im, [exp_dir 'JPEGImages/' output_img_name '.jpg']);
                    close all;
                    break;
                end                             
                                
                if 0
                    % visualize transformed mesh and where SFM keypoints fall
                    figure;
                    colormap gray;
                    axis equal;
                    hold on;
                    cmap = jet(n_kp);

                    if 0 % original
                        trisurf(nodes', xyz(1,:), xyz(2,:), xyz(3,:));
                        hold on; 
                    else % registered with Shape
                        trisurf(nodes', transf_xyz(1,:), transf_xyz(2,:), transf_xyz(3,:)); 
                        hold on;
                    end

                    for k=1:size(Shape,2)
                        plot3(Shape(1,k), Shape(2,k), Shape(3,k),  'o', 'MarkerSize', 18, 'LineWidth', 10, 'Color', cmap(k,:));
                        hold on;
                    end

                    transf_kp = bsxfun(@plus, global_s*global_R*kp3d, global_T);
                    for k=1:size(transf_kp,2)
                        plot3(transf_kp(1,k), transf_kp(2,k), transf_kp(3,k),  'o', 'MarkerSize', 18, 'LineWidth', 10, 'Color', cmap(k,:));
                        hold on;
                    end
                end

                     
                %%% find out what's visible %%%
                % orthographic projection
                P = [1 0 0 0;
                    0 1 0 0;
                    0 0 0 1];                
                rendering = RenderTriMex(P, sel_dims(2), sel_dims(1), cam_xyz, [], uint32(nodes-1))';
                
                if 0
                    %
                    % compute mask
                    %
                    % this sometimes fails when there's a very thin surface
                    % (eg. in a dining table)
                    
                    mask = rendering>0;                
                else
                    % just threshold the other mesh rendering
                    mask = render_im(:,:,2)<255;                                        
                end
                
                if 0
                    % retain largest connected component 
                    mask = getLargetCC(mask);                                
                end
                 
                masks = {mask}; % compatibility with existing stuff
                save([mask_dir output_img_name '.mat'], 'masks');
                
                % generate Pascal VOC ground truth segmentation annotations
                I_obj = uint8(mask);
                I_cls = uint8(mask);
                I_cls(mask) = uint8(h);

                cmap = VOClabelcolormap(256);
                I_obj = I_obj;
                I_cls = I_cls;
                imwrite(I_obj, cmap, [exp_dir 'SegmentationObject/' output_img_name '.png']);
                imwrite(I_cls, cmap, [exp_dir 'SegmentationClass/' output_img_name '.png']);                
                
                %
                % compute keypoint projections
                %
                transf_kp = bsxfun(@plus, global_s*global_R*kp3d, global_T);
                cam_transf_kp = bsxfun(@plus,camR*transf_kp,camT);

                %
                % check if the faces containing the keypoints are
                % visible
                
                % first assign each keypoint to closest face(s)
                face_vis = false(size(kp_obj,1),1);
                face_id = double(face_vis);
                
                nodes = nodes;
                xyz = xyz;
                %parfor k=1:size(kp_obj,1)
                                
                parfor k=1:size(kp_obj,1)             
                    if(any(isnan(kp3d(:,k))))
                        % keypoint does not exist on the mesh (due to
                        % intra-class variation)
                        face_vis(k) = false;
                        continue;
                    end
                    
                    B = cartToBary(tri_rep, (1:size(nodes,2))', repmat(kp3d(:,k),1,size(nodes,2))');                    
                    candidates = find(all((B>-eps)'));
                    
                    
                    dist = zeros(size(candidates,2),1);
                    for l=1:numel(candidates)
                        [dist(l)] = pointTriangleDistance(xyz(:,nodes(:,candidates(l)))', kp3d(:,k));
                    end
                    [val,id] = min(dist);
                    face_id(k) = candidates(id);
                    
                    if 0
                        % visualize face and point
                        trisurf(nodes(:,face_id(k))', xyz(1,:), xyz(2,:), xyz(3,:));
                        hold on;
                        plot3(kp3d(1,k), kp3d(2,k), kp3d(3,k),  'ob', 'MarkerSize', 10, 'LineWidth', 3);
                    end
                    
                    % see if face is occluded
                    face_vis(k) = any(rendering(:)==face_id(k));
                end
                
                % sometimes humans annotate points that are occluded                
                % we'll transfer the information about which keypoints were marked as visible 
                % from the pascal image 
                transf_visible_kp = sum(kp_obj,2)~=0;
                    
                visible_kp = logical(face_vis);                
                
                % there may be noise in the camera so we augment 
                % human keypoints with those that are geometrically visible
                % given the mesh and the camera
                visible_kp = visible_kp | transf_visible_kp;
                                
                if 0
                    sc(mask); hold on; 
                    %trisurf(nodes',cam_xyz(1,:),cam_xyz(2,:),cam_xyz(3,:),'FaceColor','red'); hold on;
                    plot(cam_transf_kp(1,visible_kp), cam_transf_kp(2,visible_kp),  'og', 'MarkerSize', 10, 'LineWidth', 5);
                    %plot(cam_transf_kp(1,:), cam_transf_kp(2,:), 'og', 'MarkerSize', 10, 'LineWidth', 5);
                end
                
                kp = cam_transf_kp;
                kp(:,~visible_kp) = 0;                
                                
                % save keypoints (same format as elsewhere, for compatibility)                                
                rec.objects(1).keypoints = kp'; % all 3 dims
                rec.objects(1).visib = visible_kp;                
                save([exp_dir 'merged_Correspondences_GT_BRKL/' output_img_name '.mat'], 'rec');

                % store mesh in camera frame
                tri = nodes;
                vertices = cam_xyz;
                save([cam_mesh_dir output_img_name '.mat'], 'tri', 'vertices', 'global_R', 'global_s', 'global_T', 'P', 'cam_transf_kp', 'camR', 'camT', 'visible_kp');                                 
            end        
        end  
        close all;
        time_mesh = toc(t_mesh)
    end
    
    if 1 % for running the load function
        mkdir('Dataset/Renderings/JPEGImages/');
        system('cp VOC/JPEGImages/*view* ./Dataset/Renderings/JPEGImages/');
        
        mkdir('Dataset/Renderings/MySegmentsMat/ground_truth/');
        system('cp VOC/MySegmentsMat/ground_truth/*view* ./Dataset/Renderings/MySegmentsMat/ground_truth/');
        
        mkdir('Dataset/Renderings/merged_Correspondences_GT_BRKL/');
        system('cp VOC/merged_Correspondences_GT_BRKL/*view* ./Dataset/Renderings/merged_Correspondences_GT_BRKL/');
        
        system(['cp VOC/ImageSets/Segmentation/' imgset_name '.txt ./Dataset/Renderings/']);
        
        mkdir('Dataset/Renderings/SegmentationObject/');
        system('cp VOC/SegmentationObject/*view* ./Dataset/Renderings/SegmentationObject/');

        mkdir('Dataset/Renderings/SegmentationClass/');
        system('cp VOC/SegmentationClass/*view* ./Dataset/Renderings/SegmentationClass/');  
        
        
        mkdir('Dataset/Renderings/MyMeshes/ground_truth/');
        system('cp VOC/MyMeshes/ground_truth/*view* ./Dataset/Renderings/MyMeshes/ground_truth/');
        
        d = dir('./Dataset/Renderings/MyMeshes/ground_truth/*.mat');
        
        for i=1:numel(d) % do not save the meshes here to save space
            load(['./Dataset/Renderings/MyMeshes/ground_truth/' d(i).name]); 
            save(['./Dataset/Renderings/MyMeshes/ground_truth/' d(i).name], 'global_R','global_s', 'global_T', 'P', 'cam_transf_kp', 'camR', 'camT', 'visible_kp');
        end
    end
end