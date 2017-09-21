function convert_dataset_kp_format()  
  exp_dir = add_all_paths();
  
  mask_type = 'ground_truth';
  imgset = 'all_gt_segm';

  warning off;

  sb = SegmBrowser(exp_dir, mask_type, imgset);

  dest_dir = [exp_dir 'merged_Correspondences_GT_BRKL/'];
  if(~exist(dest_dir, 'dir'))
      mkdir(dest_dir);
  end
  
  for i=1:numel(sb.categories)
      var = load(['annotations/unique_kp_' sb.categories{i} '.mat']);
      kp_names_class{i} = var.un_names;
  end

  VOCinit();
  for i=1:numel(sb.img_names)
      rec.objects = [];
      i
      if(exist([dest_dir sb.img_names{i} '.mat'], 'file'))
        %continue;
      end

      GT_rec_ours = load([exp_dir 'SegmentEval/ground_truth/overlap/' sb.img_names{i} '.mat']);
      GT_rec_ours = GT_rec_ours.Quality{1};            

      for k=1:numel(GT_rec_ours)
        rec.objects(k).duh = [];
      end

      masks = sb.get_masks(sb.img_2_whole_ids{i});
      mask_bboxes = sb.get_bboxes_from_masks(masks);
      x_min = [mask_bboxes(:).xmin];
      y_min = [mask_bboxes(:).ymin];
      x_max = [mask_bboxes(:).xmax];
      y_max = [mask_bboxes(:).ymax];
      %width = [mask_bboxes(:).xmax] - [mask_bboxes(:).xmin];
      %height = [mask_bboxes(:).ymax] - [mask_bboxes(:).ymin];
      %mask_bboxes = [x_min; y_min; width; height]'; % each row specifies a bbox
      mask_bboxes = [x_min; y_min; x_max; y_max]';

      %assert(numel([GT_rec_ours.object]) == numel([GT_rec.objects]))
      if 0
          obj_ids = 1;
          sb.get_mask_bboxes(obj_ids);
      end

      ids_checklist = false(numel(GT_rec_ours),1);

      classes_img = [sb.img_overlaps{i}.class];
      un_classes_img = unique(classes_img);
      img_name = sb.img_names{i};
      counter = 1;
      
      kp_struct = [];
      for j=1:numel(un_classes_img)
          dir_names = ['annotations/' sb.categories{un_classes_img(j)} '/' img_name '*.xml'];
          filenames = dir(dir_names);
          if(isempty(filenames))
              % can't find files for that class
              %disp(['missing files for ' sb.categories{un_classes_img(j)} ' and img ' img_name ' !']);
              continue;
          end
          
          for k=1:numel(filenames)
              kp_s = read_brkl_kp(['annotations/' sb.categories{un_classes_img(j)} '/' filenames(k).name]);
              kp_s.class = un_classes_img(j);
              kp_struct{counter} = kp_s;
              counter = counter + 1;
          end
      end
      
      LARGE_COST = 1000;
      o = LARGE_COST * ones(size(mask_bboxes,1), numel(kp_struct));      
      % match segments with keypoint objects
      for k=1:numel(kp_struct)
          kp_s = kp_struct{k};
          kp_bbox = [kp_s.visible_bounds(1) ...
              kp_s.visible_bounds(2) ...
              kp_s.visible_bounds(4) + kp_s.visible_bounds(1)...
              kp_s.visible_bounds(3) + kp_s.visible_bounds(2)];
      
          o(:,k) = -boxoverlap(mask_bboxes, kp_bbox);
          
          % wrong class get negative o, just in case      
          o([GT_rec_ours.class] ~= kp_s.class, k) = LARGE_COST;
      end
      
      
      if(size(o,2) > size(o,1))
          % more kp than masks
          o = [o; LARGE_COST * ones(size(o,2)-size(o,1), size(o,2))];
      elseif(size(o,1) > size(o,2))
          % more masks than kp
          o = [o LARGE_COST * ones(size(o,1), size(o,1)-size(o,2))];
      end
      o = o - min(o(:));
      [ass, cost] = hungarian(o);
      
      ind_costs = o(sub2ind(size(o), ass, 1:numel(ass)));
      
      if 0
          kp_s = kp_struct{1}
          kp_bbox = [kp_s.visible_bounds(2) ...
              kp_s.visible_bounds(1) ...
              kp_s.visible_bounds(2) + kp_s.visible_bounds(4) ...
              kp_s.visible_bounds(1) + kp_s.visible_bounds(3)];
          vb = kp_s.visible_bounds;
          % visualization
          sb.show_imgs(i); hold on;
          for m=1:size(mask_bboxes)
              rectangle('Position', [mask_bboxes(m,1) mask_bboxes(m,2) ...
                  mask_bboxes(m,3)-mask_bboxes(m,1) mask_bboxes(m,4)-mask_bboxes(m,2)], 'LineWidth', 5); hold on;
          end
          rectangle('Position', [vb(1) vb(2) vb(4) vb(3)], 'LineWidth', 5, 'EdgeColor', [1 0 0]); hold on;
      end
      
      for k=1:min(numel(rec.objects), numel(kp_struct))                          
          id = ass(k);
          
          % remove any element with LARGE_COST (class swaps!)
          if(ind_costs(k) >= LARGE_COST)
              continue;
          end
          
          if(isfield(kp_struct{k}, 'visib'))
              rec.objects(id).visib = kp_struct{k}.visib;
          else
              rec.objects(id).visib = [];
          end
          
          rec.objects(id).kp_names = kp_struct{k}.kp_name;
          if(isfield(kp_struct{k}, 'kp_coords'))
              n_kp = numel(kp_names_class{kp_struct{k}.class});
              
              rec.objects(id).keypoints = zeros(2,n_kp);
              [ids,ids_a,ids_b] = intersect(kp_names_class{kp_struct{k}.class}, kp_struct{k}.kp_name);
              rec.objects(id).keypoints(:,ids_a) = kp_struct{k}.kp_coords(1:2,ids_b);
              rec.objects(id).keypoints = rec.objects(id).keypoints';
              
              assert(size(rec.objects(id).keypoints ,1) == n_kp);
          else
              rec.objects(id).keypoints = [];
          end
          
          if(isfield(kp_struct{k}, 'subclass'))
              rec.objects(id).subclass =  kp_struct{k}.subclass;
          else
              rec.objects(id).subclass = [];
          end
          
          ids_checklist(id) = true;
      end
      
      for k=1:numel(ids_checklist)
          if(~ids_checklist(k))
              rec.objects(k).keypoints = [];
          end
      end
      
      save([dest_dir sb.img_names{i} '.mat'], 'rec');
  end
end

% now you can call collect...