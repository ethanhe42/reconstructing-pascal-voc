classdef CorrespBRKL
    properties(Access = private)                
        exp_dir
        corr_dir
        imgset
    end
    
    properties(Access = public)
        nPBM
                
        obj_keypoints
        img_ann % annotations for each image        
        
        warping_type
    end
    
    methods
        function obj = CorrespBRKL(exp_dir, imgset, nPBM)
            the_dir = 'merged_Correspondences_GT_BRKL/';
            obj.exp_dir = exp_dir;
            obj.corr_dir = [obj.exp_dir the_dir];
            obj.imgset = imgset;
                        
            if(~exist(obj.corr_dir, 'dir'))
                mkdir(obj.corr_dir);
            end
            
            var = load([obj.corr_dir obj.imgset '.mat']);
            obj.obj_keypoints = var.obj_keypoints;
            obj.img_ann = var.rec;
            
            assert(numel(obj.img_ann) == numel(nPBM.img_names));
            
            % assemble subcategory information
            %subcategories = cell(numel(obj.obj_keypoints),1);
            %subcategories = {};
            %for i=1:numel(nPBM.img_names)
            % subcategories = [subcategories {obj.img_ann(i).objects(:).subclass}];              
            %end
            
            obj.nPBM = nPBM;
            
            obj.warping_type = 'bbox_non_isotropic'; % 'bbox_isotropic'
        end                

        function annotate_keypoints(obj, category)
            img_names = textread([obj.exp_dir 'ImageSets/Segmentation/' obj.imgset '.txt'], '%s');

%             if 0
%                 mkdir([exp_dir '/JPEGImages/']);
%                 for i=1:numel(img_names)
%                     copyfile([obj.exp_dir 'JPEGImages/' img_names{i} '.jpg'], ['./VOC2010/JPEGImages/']);
%                 end
%             end

            % show all training examples
            if 0
                obj.nPBM.show_all_category_imgs(obj.exp_dir, category, img_names)
            end

            % retrieve all examples
            VOCInit();
            category_corr = [];
            for i=1:numel(img_names)
                close all;
                img_names{i}
                filename = [obj.corr_dir img_names{i} '.mat'];

                I = imread([obj.exp_dir 'JPEGImages/' img_names{i} '.jpg']);

                rec=PASreadrecord(sprintf(VOCopts.annopath,img_names{i}));
                change_count = 0;
                for j=1:numel(rec.objects)
                    o = rec.objects(j);
                    if(strcmp(o.class, category))
                        if(~exist(filename))
                            c = obj.annotate_keypoints_in_bbox(I,rec.objects(j).bndbox, category);
                            rec.objects(j).keypoints = c;
                            change_count = change_count + 1;
                        elseif 0 % show current annotations
                            load(filename);
                            obj.show_corresp_category_bbox(I, rec.objects(j).keypoints, rec.objects(j).bndbox, category);
                            pause;
                        end
                    end
                end
                
                if(change_count>0)
                    save(filename, 'rec');
                end
            end          
        end

        function c = annotate_keypoints_in_bbox(obj, I, bndbox, category)
            [keypoint_names, colors] = obj.get_category_keypoint_names(category);

            range_x = bndbox.xmin:bndbox.xmax;
            range_y = bndbox.ymin:bndbox.ymax;
            I = I(range_y, range_x,:);

            imshow(I);
            for i=1:numel(keypoint_names)
                keypoint_names{i}
                [p_x, p_y, key] = ginput(1)

                if(key == 1) % left mouse button
                    hold on;
                    plot(p_x, p_y, 'x', 'MarkerSize', 10, 'LineWidth', 3, 'MarkerEdgeColor', colors(i,:));
                    p_x = p_x + range_x(1);
                    p_y = p_y + range_y(1);

                    kp{i} = [p_x p_y];
                else
                    kp{i} = [0 0];
                end
            end
            hold off;

            kp
            c = cell2mat(kp');
        end

        function [keypoint_names,colors] = get_category_keypoint_names(obj, category)
            if(isnumeric(category))
                category = obj.nPBM.categories{category};
            end
            
            var = load([obj.exp_dir 'Browser/unique_kp_' category '.mat']);
            keypoint_names = var.un_names;
            %colors = lines(numel(keypoint_names));
            colors = obj.get_cmap_and_symbols();
        end

        function [kps, kps_obj_ids] = get_keypoints_in_masks(obj, mask_ids, whole_or_parts, all_or_dominant_obj)
            DefaultVal('*whole_or_parts', '''whole''');
            DefaultVal('*all_or_dominant_obj', '''dominant_obj''');

            kps = [];

            % find the image the mask is from, and which objects are there
            if(strcmp(whole_or_parts, 'whole'))
                img_ids = obj.nPBM.whole_2_img_ids(mask_ids);
            else
                whole_ids = obj.nPBM.part_2_whole_ids(mask_ids);                
                img_ids = obj.nPBM.whole_2_img_ids(whole_ids);
            end

            %un_img_ids = unique(img_ids);
            un_img_ids = img_ids;
            obj_ids = cell(numel(un_img_ids),1);
            for i=1:numel(img_ids)
                obj_ids{i} = find(obj.nPBM.obj_2_img_ids==un_img_ids(i));
            end

            masks = obj.nPBM.get_masks(mask_ids, whole_or_parts);
            kps = cell(numel(masks),1);
            kps_obj_ids = kps;
            for i=1:numel(masks)       
                % grow the mask a few pixels
                masks{i} = imdilate(masks{i}, ones(4,4));
                           
                obj_counter = 1;   
                n_present = [];
                for j=1:numel(obj_ids{i})
                    if(isempty(obj.obj_keypoints{obj_ids{i}(j)}))
                        continue;
                    end
                    
                    present = (obj.obj_keypoints{obj_ids{i}(j)}(:,1)~=0);
                    n_present(obj_counter) = sum(present);
                    if(n_present(obj_counter)==0)                        
                        continue;
                    end
                    
                    in_masks = zeros(size(present));
                    
                    % account foir noisy keypoints set outside the image
                    these_kps = round(obj.obj_keypoints{obj_ids{i}(j)}(present,:));
                    these_kps(these_kps(:,1)>size(masks{i},2), 1) = size(masks{i},2);
                    these_kps(these_kps(:,2)>size(masks{i},1), 2) = size(masks{i},1);
                    these_kps(1>these_kps) = 1;
                    in_masks(present) = masks{i}(sub2ind(size(masks{i}), these_kps(:,2), these_kps(:,1)));
                    
                    if 0 % debug
                        %sc(masks{i});
                        obj.nPBM.show_imgs(img_ids(i));
                        hold on;
                        scatter(these_kps(:,1), these_kps(:,2));
                    end
                    
                    if(sum(in_masks) == 0)
                        n_present(obj_counter) = [];
                        continue;
                    end
                    kps{i}{obj_counter} = obj.obj_keypoints{obj_ids{i}(j)};
                    kps{i}{obj_counter}((~in_masks),:) = nan;
                    kps_obj_ids{i}(obj_counter) = obj_ids{i}(j);
                    obj_counter = obj_counter + 1;                    
                end
                
                if(strcmp(all_or_dominant_obj, 'dominant_obj'))
                    if(isempty(n_present))
                        continue;
                    end
                    [m, id] = max(n_present);
                    kps{i} = kps{i}{id};
                    kps_obj_ids{i} = kps_obj_ids{i}(id);
                    % only retain the keypoints from the dominant object
                end                
            end
        end

        function [aspect_ids, distortion] = cluster_aspects_D(obj, n_aspects, D)
            aspect_ids = [];
            aspect_prototypes = [];
            distortion = [];

            if 0 %%%% Debug %%%%                                
                % visualize sorted 5 nearest neighbors for 3 examples
                %rp = randperm(size(D,1));
                %D_dbg = D(rp(1:3),:);
                D_dbg = D([10 20 30], :);
                for i=1:size(D_dbg,1)
                    [dist,closest_n] = sort(D_dbg(i,:), 'ascend');
                    closest_n = closest_n(1:5);
                    dist_n = dist(1:5);
                    for j=1:numel(dist_n)
                        titles{j} = sprintf('%f', dist_n(j));
                    end

                    obj.show_obj_imgs_and_keypoints(obj_ids(closest_n), titles);
                end
            end

            if 1 % kmedoids
                [aspect_prototypes, distortion, aspect_ids] = kmedoid(D, n_aspects, 100000);
            elseif 0 % agglomerative clustering
                % ward seems the least bad, none seems good for low
                % n_aspects
                Z = linkage(squareform(D,'tovector'), 'ward');
                % dendrogram(Z)
                
                aspect_ids = cluster(Z, 'MaxClust', n_aspects);
                hist(aspect_ids)
                distortion = [];
            elseif 1 % spectral
                a  = Dgraph(exp(-0.45*D));
                bin_mat = a.spectral_cuts('yu', 'NSegments', n_aspects);
                [duh,aspect_ids] = find(bin_mat);
            end
        end

        function [aspect_ids, distortion] = cluster_aspects(obj, n_aspects, coords, obj_ids)
            aspect_ids = [];
            aspect_prototypes = [];
            distortion = [];

            D = distance_corresp_lists(coords, coords);

            if 0 %%%% Debug %%%%                                
                % visualize sorted 5 nearest neighbors for 3 examples
                %rp = randperm(size(D,1));
                %D_dbg = D(rp(1:3),:);
                D_dbg = D([10 20 30], :);
                for i=1:size(D_dbg,1)
                    [dist,closest_n] = sort(D_dbg(i,:), 'ascend');
                    closest_n = closest_n(1:5);
                    dist_n = dist(1:5);
                    for j=1:numel(dist_n)
                        titles{j} = sprintf('%f', dist_n(j));
                    end

                    obj.show_obj_imgs_and_keypoints(obj_ids(closest_n), titles);
                end
            end

            if 1 % kmedoids
                [aspect_prototypes, distortion, aspect_ids] = kmedoid(D, n_aspects, 5000000);
            elseif 0 % agglomerative clustering
                % ward seems the least bad, none seems good for low
                % n_aspects
                D = D - diag(diag(D));
                Z = linkage(squareform(D,'tovector'), 'ward');
                %Z = linkage(squareform(D,'tovector'));
                % dendrogram(Z)
                
                aspect_ids = cluster(Z, 'MaxClust', n_aspects);
                hist(aspect_ids)
                distortion = [];
            elseif 0 % spectral
                a  = Dgraph(exp(-0.45*D));
                bin_mat = a.spectral_cuts('yu', 'NSegments', n_aspects);
                [duh,aspect_ids] = find(bin_mat);
            end
        end

        function [f_coords, frame_ids] = get_ref_frames(obj, obj_kp, obj_mask)
          obj_kp_bbox = obj.img2bbox_coords(obj_kp, obj_mask);
          c = obj_kp_bbox{1};
          def = find(c(:,1)~=0);
          
          c = c(def,:)';
          %obj_D = pdist2(c',c','euclidean');
          %obj_D = obj_D + 10000*eye(size(obj_D));
          %[d, id] = min(obj_D);
          %all_ref_frames = [(1:numel(d))' id'];
          [a,b] = find(triu(ones(size(c,2), size(c,2)),1));
          all_ref_frames = [a b];
          
          % reference frame is defined by two points
          f_coords = [c(:,all_ref_frames(:,1)); c(:,all_ref_frames(:,2))];

          frame_ids = def(all_ref_frames);
          if 0 
            obj.plot_ref_frames(f_coords);
          end
        end
        
        function [proj_coords, T] = proj_ref_frames(obj, ref_f_coords, p, T)
          DefaultVal('*T', '[]');
          
          if(isempty(T))
            % find transformation between bbox and new coord frames
            base_points = [0 0; 0 1];
            T = cell(size(ref_f_coords,2),1);
            for i=1:size(ref_f_coords,2)            
              new_points = reshape(ref_f_coords(:,i), 2,2)';
              T{i} = cp2tform(base_points, new_points, 'nonreflective similarity');
            end
          end
          
          proj_coords = cell(numel(T),1);
          for i=1:numel(T)
            proj_coords{i} = tforminv(T{i}, p');
          end
        end
        
        function [proj_coords, T] = proj_ref_frames_fast(obj, ref_f_coords, p, T)
          DefaultVal('*T', '[]');
          
          if(isempty(T))
            % find transformation between bbox and new coord frames
            base_points = [0 0; 0 1];
            T = cell(size(ref_f_coords,2),1);
            for i=1:size(ref_f_coords,2)
              T{i} = compute_T_nonref_sim(ref_f_coords(:,i), base_points, 'inv');
            end
          end
          
          proj_coords = cell(numel(T),1);
          for i=1:numel(T)
            X = p';
            M = T{i};
            X1 = [X ones(size(X,1),1)];   % Convert X to homogeneous coordinates
            U1 = X1 * M;                  % Transform in homogeneous coordinates
            proj_coords{i}  = U1(:,1:end-1); % Convert homogeneous coordinates to U   
          end
        end
        
        function [inv_proj_coords, T] = inv_proj_ref_frames(obj, ref_f_coords, proj_p, T)
          DefaultVal('*T', '[]');
          
          if(isempty(T))
            % find transformation between bbox and new coord frames
            base_points = [0 0; 0 1];
            p = zeros(2, size(ref_f_coords,2));
            T = cell(size(ref_f_coords,2),1);
            for i=1:size(ref_f_coords,2)
              new_points = reshape(ref_f_coords(:,i), 2,2)';
              T{i} = cp2tform(base_points, new_points, 'nonreflective similarity');
            end
          end
          
          inv_proj_coords = cell(numel(T),1);
          for i=1:numel(T)
            inv_proj_coords{i} = tformfwd(T{i}, proj_p');
          end
        end
        
        function [inv_proj_coords, T] = inv_proj_ref_frames_fast(obj, ref_f_coords, proj_p, T)
          DefaultVal('*T', '[]');
          
          if(isempty(T))
            % find transformation between bbox and new coord frames
            base_points = [0 0; 0 1];
            p = zeros(2, size(ref_f_coords,2));
            T = cell(size(ref_f_coords,2),1);
            for i=1:size(ref_f_coords,2)
              T{i} = compute_T_nonref_sim(ref_f_coords(:,i), base_points, 'direct');
            end
          end
          
          inv_proj_coords = cell(numel(T),1);
          for i=1:numel(T)
            X = proj_p';
            M = T{i};
            X1 = [X ones(size(X,1),1)];   % Convert X to homogeneous coordinates
            U1 = X1 * M;                  % Transform in homogeneous coordinates
            inv_proj_coords{i}  = U1(:,1:end-1); % Convert homogeneous coordinates to U         
          end
        end
        
        function plot_ref_frames(obj, ref_frames)
          error('fix');
            hold on;
            
            for i=1:size(ref_frames,2)
                plot(ref_frames(1,i), ref_frames(2,i), 'ob');
                arrow(ref_frames(1:2,i), ref_frames(1:2,i)+ref_frames(3:4,i));
                arrow(ref_frames(1:2,i), ref_frames(1:2,i)+ref_frames(5:6,i));
            end
        end
        
        %%% coordinate normalization %%%
        function [n_coords, warp_pars] = warp_coords(obj, coords, masks)
            n_coords = cell(1,numel(coords));
            warp_pars = n_coords;
            sz = cellfun(@(a) size(a,1), coords);
            
            coords_present = false(max(sz), numel(coords));
            
            %switch(obj.warping_type)
            %    case {'bbox_isotropic', 'bbox_non_isotropic'}
            %        bboxes = obj.nPBM.get_bboxes(obj_ids);
            %end
                                    
            for i=1:numel(coords)  
                if(isempty(coords{i}))
                    n_coords{i} = [];
                    warp_pars{i} = [];
                    continue;
                end
                
                coords_present(:,i) = (coords{i}(:,1)~=0);
                                
                n_coords{i} = coords{i};
                missing = ~coords_present(:,i);                
                n_coords{i}(missing,:) = nan(sum(missing),2);
                
                switch(obj.warping_type)                                     
                    case {'bbox_isotropic', 'bbox_non_isotropic'}                                              
                        [ids_y, ids_x] = find(masks(:,:,i));
                        max_x = max(ids_x); min_x = min(ids_x);
                        max_y = max(ids_y); min_y = min(ids_y);
                        center_bbox_x = mean([min_x max_x]);
                        center_bbox_y = mean([min_y max_y]);
                        width = max(1, max_x - min_x);
                        height = max(1,max_y - min_y);
                            
                        %bndboxes = obj.nPBM.get_bboxes_from_masks(masks);
                        
                        if(strcmp(obj.warping_type, 'bbox_isotropic'))
                            scaling = max(height, width);
                            n_coords{i}(coords_present(:,i),2) = -(coords{i}(coords_present(:,i),2) - center_bbox_y)/scaling; % want cartesian frame to grow up not down
                            n_coords{i}(coords_present(:,i),1) = (coords{i}(coords_present(:,i),1) - center_bbox_x)/scaling;
                            warp_pars{i} = scaling;
                        elseif(strcmp(obj.warping_type, 'bbox_non_isotropic'))
                            warp_pars{i} = [width height];
                            n_coords{i}(coords_present(:,i),2) = -(coords{i}(coords_present(:,i),2) - center_bbox_y)/height; % want cartesian frame to grow up not down
                            n_coords{i}(coords_present(:,i),1) = (coords{i}(coords_present(:,i),1) - center_bbox_x)/width;
                            if 0 % debug
                                [Imgs, tl_bb_corners] = obj.nPBM.get_object_Imgs(obj_ids(i), true);
                                sc(Imgs{1});
                            end
                        end
                    case{'none'}
                        n_coords{i} = coords{i};
                end
                                  
                if 0 % debug    
                    axis xy;
                    axis([-2 2 -2 2]);
                    hold on;
                    plot(n_coords{i}(coords_present(:,i),1), n_coords{i}(coords_present(:,i),2), 'og');
                    %sc(masks{i}); hold on;                        
                    %plot(coords{i}(coords_present(:,i),1), coords{i}(coords_present(:,i),2), 'ob');
                end
            end
        end
        
        function [bbox_coords, bbox] = img2bbox_coords(obj, coords, masks)
            for i=1:numel(masks)
                assert(numel(coords) == numel(masks));
                [ids_y, ids_x] = find(masks{i});
                max_x = max(ids_x); min_x = min(ids_x);
                max_y = max(ids_y); min_y = min(ids_y);
                
                no_change = coords{i}(:,1)==0;
                
                bbox = [min_x min_y max_x max_y];
                bbox_coords{i} = coords{i};
                bbox_coords{i}(~no_change,:) = bbox_coords{i}(~no_change,:) - repmat([min_x min_y], sum(~no_change), 1);
            end
        end
        
        %%%%%%%%%%%%%
        % Visualization  %
        %%%%%%%%%%%%%

        function show_all_inst_of_keypoint(obj, keypoint_ids, coords, obj_ids, names, line_segment_kp)            
            % show keypoints from multiple examples overlaid
            DefaultVal('*names', '[]');
            DefaultVal('*with_lines', 'true');

            figure;
            %axis manual;
            axis([-0.6 .6 -0.6 0.6]);
            set(gca, 'PlotBoxAspectRatio', [1 1 1])
            hold on;                    

            [cmap, symbols] = obj.get_cmap_and_symbols();
            for i=1:numel(coords)
                obj.plot_keypoints(keypoint_ids, coords{i}, symbols, cmap);
            end

            if(~isempty(names))          
                [a,b] = find(~isnan(cell2mat(coords)));
                [keypoints_present, pos] = unique(a);
                
                [legend_h,object_h,plot_h,text_strings] = legend(names(keypoints_present), 'Interpreter', 'none') 
                
                for j=1:numel(keypoints_present)
                    set(plot_h(j), 'MarkerEdgeColor', cmap(keypoint_ids(keypoints_present(j)),:), 'Marker', symbols(keypoint_ids(keypoints_present(j))));
                end            
                
                % it's stupid i know, but set above changes a small number
                % of symbols...
                duh = get(gca);
                delete(duh.Children)
                for i=1:numel(coords)
                    obj.plot_keypoints(keypoint_ids, coords{i}, symbols, cmap);
                end
            end

            for i=1:numel(coords)
                if(isempty(coords{i}))
                    continue;
                end
                
                present = (coords{i}(line_segment_kp,1)~=0);
                ids_present = find(present);
                c = coords{i}(line_segment_kp(present),:);
                
                %%% compute minimum spanning tree between selected keypoints
                %D = pdist2(c,c, 'euclidean');
                %m = mst(sparse(D));
                %[e_i, e_j] = find(triu(m));
                
                if(numel(ids_present)>1 && with_lines)
                    [e_i, e_j] = cartprod_mex((1:numel(ids_present))', (1:numel(ids_present))')
                    jl_plot_lines([c(e_i(:),2) c(e_i(:), 1) c(e_j(:), 2) c(e_j(:),1)]');
                end
            end

            axis xy;

            handle = uicontrol(gcf, 'Style', 'pushbutton', 'String', 'Show Keypoint Image');
            set(handle, 'Position', [210 20 200 50]);            

            set(handle, 'Callback', @(src, event) ginput_wrapper(obj, src, event, obj_ids, coords));
        end

        function plot_keypoints(obj, keypoint_ids, coords, symbols, cmap)
            % plot keypoints for a single example
            hold on;
            for j=1:numel(keypoint_ids)                
                if(~isempty(coords) && coords(keypoint_ids(j),1)~=0)
                    plot(coords(keypoint_ids(j),1), coords(keypoint_ids(j),2), symbols(keypoint_ids(j)), 'LineWidth', 10, 'MarkerSize', 10, 'MarkerEdgeColor', cmap(keypoint_ids(j),:));
                end
            end
        end
                
        function show_obj_imgs_and_keypoints(obj, obj_ids, titles)
            DefaultVal('*titles', '[]');
            %obj.obj_keypoints
            
            coords = obj.obj_keypoints(obj_ids);
            
            just_bbox = true;
            [Imgs, top_left_bbox_corners] = obj.nPBM.get_object_Imgs(obj_ids, just_bbox);
            % load all images
            
            handles = montage_new(Imgs, titles);
            
            [cmap, symbols] = obj.get_cmap_and_symbols();
                        
            for i=1:numel(handles)  
                if(isempty(coords{i}))
                    continue;
                end
                
                present = (coords{i}(:,1)~=0);
                
                axes(handles(i));
                coords{i}(present,1) = coords{i}(present,1) - top_left_bbox_corners(1,i);
                coords{i}(present,2) = coords{i}(present,2) - top_left_bbox_corners(2,i);
                obj.plot_keypoints(1:size(coords{i},1), coords{i}, symbols, cmap);
            end            
        end
        
        function plot_keypoints_obj(obj, obj_id)
            coords = obj.obj_keypoints{obj_id};
           
            [cmap, symbols] = obj.get_cmap_and_symbols();
            
            present = (coords(:,1)~=0);
            
            obj.plot_keypoints(1:size(coords,1), coords, symbols, cmap);
        end
        
        function handles = show_wholes_and_keypoints(obj, whole_ids, show_wholes)
            DefaultVal('*show_wholes', 'true');
            show_kp_ids = true;
            %[kps, kps_obj_ids] = obj.get_keypoints_in_masks(whole_ids);
            kps = obj.obj_keypoints(whole_ids);
            
            %if(any(cellfun(@isempty, kps)))                
            %    handles = -1;
            %    return;
            %end
            
            Imgs = obj.nPBM.get_Imgs(obj.nPBM.whole_2_img_ids(whole_ids));
            masks = obj.nPBM.get_masks(whole_ids);            
            
            bboxes = obj.nPBM.get_bboxes_from_masks(masks);
                     
            top_left_bbox_corners = zeros(2, numel(bboxes));
            for i=1:numel(bboxes)
                Imgs{i} = Imgs{i}(bboxes(i).ymin:bboxes(i).ymax,bboxes(i).xmin:bboxes(i).xmax, :);
                masks{i} = masks{i}(bboxes(i).ymin:bboxes(i).ymax, bboxes(i).xmin:bboxes(i).xmax);
                top_left_bbox_corners(:,i) = [bboxes(i).xmin; bboxes(i).ymin];
            end
            
            if(show_wholes)
                [Imgs, handles] = subplot_auto_transparent_multiple_imgs(masks, Imgs);
            else
                for i=1:numel(masks)
                    masks{i} = false;
                end
                [Imgs, handles] = subplot_auto_transparent_multiple_imgs(masks, Imgs);
            end
            
            [cmap, symbols] = obj.get_cmap_and_symbols();            
            for i=1:numel(handles)
                if(isempty(kps{i}))
                    continue;
                end
                
                present = (kps{i}(:,1)~=0);
                
                axes(handles(i));
                kps{i}(present,1) = kps{i}(present,1) - top_left_bbox_corners(1,i);
                kps{i}(present,2) = kps{i}(present,2) - top_left_bbox_corners(2,i);
                obj.plot_keypoints(1:size(kps{i},1), kps{i}, symbols, cmap);
                
                if(show_kp_ids)
                  p = kps{i}(present,:);
                  pres = find(present);
                  for j=1:size(p,1)
                    label = int2str(pres(j));
                    text(p(j,1), p(j,2), label, 'Fontsize', 15, 'Color', [1 1 1]);
                    %plot(kp_coords{i}(present,1), kp_coords{i}(present,2), 
                  end
                end
            end
        end
        
        function handles = show_masks_and_keypoints(obj, masks, kp_coords, obj_ids, local_feats, show_kp_ids)
            DefaultVal('*local_feats', '[]');
            DefaultVal('*show_kp_ids', 'false');
            [Imgs] = obj.nPBM.get_object_Imgs(obj_ids, false);
            
            bboxes = obj.nPBM.get_bboxes_from_masks(masks);
            top_left_bbox_corners = zeros(2, numel(bboxes));
            for i=1:numel(bboxes)
                Imgs{i} = Imgs{i}(bboxes(i).ymin:bboxes(i).ymax,bboxes(i).xmin:bboxes(i).xmax, :);
                masks{i} = masks{i}(bboxes(i).ymin:bboxes(i).ymax, bboxes(i).xmin:bboxes(i).xmax);
                top_left_bbox_corners(:,i) = [bboxes(i).xmin; bboxes(i).ymin];
            end
            
            [Imgs, handles] = subplot_auto_transparent_multiple_imgs(masks, Imgs);
            
            [cmap, symbols] = obj.get_cmap_and_symbols();            
            for i=1:numel(handles)
                present = (kp_coords{i}(:,1)~=0);
                
                axes(handles(i));
                kp_coords{i}(present,1) = kp_coords{i}(present,1) - top_left_bbox_corners(1,i);
                kp_coords{i}(present,2) = kp_coords{i}(present,2) - top_left_bbox_corners(2,i);
                obj.plot_keypoints(1:size(kp_coords{i},1), kp_coords{i}, symbols, cmap);
                
                if(show_kp_ids)
                  p = kp_coords{i}(present,:);
                  pres = find(present);
                  for j=1:size(p,1)
                    label = int2str(pres(j));
                    text(p(j,1), p(j,2), label, 'Fontsize', 15, 'Color', [1 1 1]);
                    %plot(kp_coords{i}(present,1), kp_coords{i}(present,2), 
                  end
                end
            end
            
            % local feats
            if(~isempty(local_feats))
                for i=1:numel(handles)
                    local_feats{i} = local_feats{i}';
                    axes(handles(i));
                    local_feats{i}(:,1) = local_feats{i}(:,1) - top_left_bbox_corners(1,i);
                    local_feats{i}(:,2) = local_feats{i}(:,2) - top_left_bbox_corners(2,i);
                    plot(local_feats{i}(:,1), local_feats{i}(:,2), 'or');
                end
            end
        end
        
        function show_corresp_category_bbox(obj, I, corr, bndbox, category)
            % get only region of interest
            range_x = bndbox.xmin:bndbox.xmax;
            range_y = bndbox.ymin:bndbox.ymax;
            I = I(range_y, range_x,:);
            
            corr
            imshow(I);
            [keypoint_names, colors] = obj.get_category_keypoint_names(category);
            for i=1:numel(keypoint_names)
                p_x = corr(i,1) - range_x(1);
                p_y = corr(i,2) - range_y(1);
                if(p_x==0) % occluded
                    continue;
                end
                hold on;
                plot(p_x, p_y, 'x', 'MarkerSize', 10, 'LineWidth', 3, 'MarkerEdgeColor', colors(i,:));
            end
            hold off;
        end 
        
        function show_categ_objects(obj, category, ids_within_categ)            
            y = obj.nPBM.get_overlaps_wholes(1:numel(obj.nPBM.whole_2_img_ids));
            [q, class] = max(y, [], 2);
            in_class = find(class == category);
            obj.show_wholes_and_keypoints(in_class(ids_within_categ))
        end
        
        function w_I = warp_img(obj, I, warp_pars)
            STD_SIZE = 200;
            switch obj.warping_type                
                case 'bbox_isotropic'
                    w_I = imresize(I, STD_SIZE/max(size(I,1), size(I,2)));
                case 'bbox_non_isotropic'
                    w_I = imresize(I, [STD_SIZE STD_SIZE]);
            end
        end
        
        function show_warped_objects(obj, obj_ids, warp_pars)
            bbox_imgs = obj.nPBM.get_object_Imgs(obj_ids, true);
            
            for i=1:numel(bbox_imgs);
                w_I{i} =warp_img(obj, bbox_imgs{i}, warp_pars{i})
            end
            
             montage_new(w_I);
        end     
        
        function show_matching_color_code(obj, obj_1, obj_2, P1_bbox, P2_bbox, show_arrows)
          DefaultVal('*show_arrows', 'false');
          P1_bbox = round(P1_bbox);
          P2_bbox = round(P2_bbox);
          
          % select bbox of best mask for each object
          Imgs =  obj.nPBM.get_Imgs(obj.nPBM.img_names(obj.nPBM.obj_2_img_ids([obj_1 obj_2])));
          [best_wholes_ids, q] = obj.nPBM.get_best_wholes([obj_1 obj_2]);
          best_masks = obj.nPBM.get_masks(best_wholes_ids, 'whole');
          bboxes = obj.nPBM.get_bboxes_from_masks(best_masks);
          I1 = Imgs{1}(bboxes(1).ymin:bboxes(1).ymax, bboxes(1).xmin:bboxes(1).xmax,:);
          I2 = Imgs{2}(bboxes(2).ymin:bboxes(2).ymax, bboxes(2).xmin:bboxes(2).xmax,:);

          %[Imgs, tl_bb_corners] = obj.nPBM.get_object_Imgs([obj_1 obj_2], true)
          %I1 = Imgs{1};
          %I2 = Imgs{2};
          HORIZ_SPACE = 50;
          
          I = zeros(max(size(I1,1),size(I2,1)), size(I1,1)+size(I2,1)+HORIZ_SPACE,3);
          
          I(1:size(I1,1), 1:size(I1,2),:) = I1;
          I(1:size(I2,1), size(I1,2)+HORIZ_SPACE:size(I1,2)+size(I2,2)+HORIZ_SPACE-1,:) = I2;

          % color varies vertically
          [p1_j,ids] = sort(P1_bbox(2,:), 'ascend');         
          
          sc(I); hold on;          
          cmap = colormap('hsv'); % 'jet', 'hsv', 'hot', 'cool'
          n_colors = size(cmap,1);
          N_PER_COLOR = max(10, ceil(numel(ids)/n_colors));
          counter = 1;
          for i=1:N_PER_COLOR:numel(ids)
            range = i:min(numel(ids),i+N_PER_COLOR-1);
            c = cmap(counter,:);
            plot(P1_bbox(1,ids(range)), P1_bbox(2,ids(range)), 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', c);
            plot(P2_bbox(1,ids(range))+size(I1,2)+1+HORIZ_SPACE, P2_bbox(2,ids(range)), 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor', c); 
            counter = counter + 1;
          end
          
          %%%% to visualize individual correspondences
          if show_arrows
            for i=1:numel(ids)
              line([P1_bbox(1,i) P2_bbox(1,i)+size(I1,2)+1+HORIZ_SPACE], [P1_bbox(2,i) P2_bbox(2,i)]);
            end
          end  
        end
    end

    methods (Access = private)
        function p = ginput_wrapper(obj, hObj, event, obj_ids, coords)
            assert(numel(coords) == numel(obj_ids));
            p = ginput(1);                

            m = zeros(numel(coords),1);
            for i=1:numel(coords)
                if(~isempty(coords{i}))
                    d = pdist2(coords{i}, p, 'euclidean');                    
                    m(i) = min(d);
                else
                    m(i) = inf;
                end
            end

            [real_m, id] = min(m);
            img_name = obj.nPBM.img_names{obj.nPBM.obj_2_img_ids(obj_ids(id))};
            I = imread([obj.exp_dir 'JPEGImages/' img_name '.jpg']);
            figure;            

            [cmap, symbols] = get_cmap_and_symbols(obj);

            keypoint_ids = find(coords{id}(:,1) ~= 0);

            obj_kp = obj.obj_keypoints{obj_ids(id)};
            subplot(1,2,1);   imshow(I); hold on;
            subplot(1,2,1); obj.plot_keypoints(keypoint_ids, obj_kp, symbols, cmap); title('original image');
            subplot(1,2,2); obj.plot_keypoints(keypoint_ids, coords{id}, symbols, cmap);  title('warped keypoints');
            set(gca, 'PlotBoxAspectRatio', [1 1 1])
            axis([-0.60 .60 -0.60 0.60]);
        end

        function [cmap, symbols] = get_cmap_and_symbols(obj)
            cmap = [     0.2260    0.2529    0.3463
                0.3610    0.8840    0.3644
                0.3246    0.1963    0.1715
                0.0836    0.1214    0.7954
                0.5127    0.5437    0.4927
                0.8329    0.3146    0.3546
                0.9046    0.3820    0.7751
                0.7236    0.7915    0.2368
                0.3830    0.8392    0.8448
                0.2980    0.6802    0.8165
                0.6917    0.4169    0.8462
                0.8805    0.6429    0.3702
                0.9245    0.2141    0.3832
                0.0813    0.6173    0.8613
                0.4827    0.6752    0.4639
                0.1283    0.6010    0.5705
                0.1283    0         0
                0.1283    0.6010    0
                0.1283    0         0.43
                0.1283    0.98      0.15
                0         0.6010    0.5705
                0         0.6010    0
                0.5       0.3       0.9];

            % symbols for plotting 
            symbols = 'oooodddd++++ssssxxxx....';
        end
    end
 end