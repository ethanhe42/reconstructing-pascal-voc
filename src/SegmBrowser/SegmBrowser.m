% September 2012, Copyright Joao Carreira 
classdef SegmBrowser    
    properties     
        VOCopts 
        
        exp_dir
        imgset
        mask_type
        
        categories                
        
        img_names        
        
        % files        
        all_whole_ids_file
        
        % hash table mapping img names to image ids
        img_name_to_id_hash
        
        % object_id 2 img ids
        obj_2_img_ids
        global_2_local_obj_ids
        
        % img 2 whole global table / whole 2 part global table
        whole_2_part_ids
        img_2_whole_ids        
        
        % tables mapping part ids to whole and whole to img ids
        part_2_whole_ids
        whole_2_img_ids
        
        % tables mapping global to local ids (e.g. id of part in whole)        
        global_2_local_whole_ids
        global_2_local_part_ids                
                
        % overlap between segments and object gt  
        img_overlaps  
    end
    
    properties (SetAccess = private)
        overlap_type
    end
    
    
    methods
        function obj = SegmBrowser(exp_dir, mask_type, imgset, overlap_type)            
            DefaultVal('*overlap_type', '''overlap''')
            obj.exp_dir = exp_dir;
            obj.imgset = imgset;
            obj.mask_type = mask_type;                        
            obj.overlap_type = overlap_type;            
                                   
            obj.img_names = textread([obj.exp_dir 'ImageSets/Segmentation/' obj.imgset '.txt'], '%s');
                        
            % right now it's formatted for VOC
            VOCinit();
            obj.VOCopts = VOCopts;
            
            for i=1:20
                obj.categories{i} = VOC09_id_to_classname(i);
            end           
            
            browser_dir = [obj.exp_dir 'Browser/' obj.mask_type '/'];
            if(~exist(browser_dir))
                mkdir(browser_dir);
            end            
            
            if(strcmp(obj.overlap_type, 'overlap'))
                obj.all_whole_ids_file = [browser_dir 'all_wholes_' imgset '.mat'];
            else
                obj.all_whole_ids_file = [browser_dir 'all_wholes_' imgset '_' obj.overlap_type '.mat'];
            end
            
            no_GT = false; % assume there are ground truth masks
            
            %
            % compute all whole ids
            %
            obj.whole_2_img_ids = [];
            global_2_local_whole_ids = [];
            if(~exist(obj.all_whole_ids_file, 'file'))
              disp('Computing and caching dataset meta-data (overlap with GT, etc)...');
              img_name_to_id_hash = zeros(numel(obj.img_names),1);
              for i=1:numel(obj.img_names)
                img_name_to_id_hash(i) = str2double(obj.img_names{i}([1:4 6:end]));
              end
            
                previous_id = 0;
                img_2_whole_ids = cell(numel(obj.img_names),1);
                whole_2_img_ids = img_2_whole_ids;
                img_overlaps = img_2_whole_ids;
                                
                for i=1:numel(obj.img_names)
                  if(exist([exp_dir 'SegmentEval/' mask_type '/' obj.overlap_type '/' obj.img_names{i} '.mat'], 'file'))
                    var = load([exp_dir 'SegmentEval/' mask_type '/' obj.overlap_type '/' obj.img_names{i} '.mat']);
                  else          
                      if(~exist([exp_dir 'SegmentEval/' mask_type '/' obj.overlap_type '/'], 'dir'))
                          mkdir([exp_dir 'SegmentEval/' mask_type '/' obj.overlap_type '/']);
                      end
                      
                      % this "try catch" command is not great. The idea is
                      % that the test set won't have ground truth and this
                      % is detected automatically. But if there are
                      % problems, the quality files will have to be removed
                      % manually.
                      try
                          parfor k=i:numel(obj.img_names)                          
                          %for k=i:numel(obj.img_names)                          
                              if(strcmp(obj.overlap_type, 'overlap_and_bbox_edges'))
                                  SvmSegm_segment_quality_all_JL_bbox(obj.img_names(k), obj.exp_dir, obj.mask_type, obj.overlap_type);
                              else
                                  SvmSegm_segment_quality_all_JL(obj.img_names(k), obj.exp_dir, obj.mask_type, obj.overlap_type);
                              end
                          end
                          var = load([exp_dir 'SegmentEval/' mask_type '/' obj.overlap_type '/' obj.img_names{i} '.mat']);
                      catch
                          % generate fake "Quality" structures
                          disp('Didn''t find GT annotations. Generating fake quality structures.');
                          Quality = [];
                          no_GT = true;
                          var = load([exp_dir 'MySegmentsMat/' obj.mask_type '/' obj.img_names{i} '.mat']);                                  

                          Quality(1).q = -ones(size(var.masks,3),1);
                          Quality(1).image_name = obj.img_names{i};
                          save([exp_dir 'SegmentEval/' mask_type '/' obj.overlap_type '/' obj.img_names{i} '.mat'], 'Quality');
                          
                          var.masks = [];
                          var.Quality = Quality;
                      end
                  end
                  
                  Quality = var.Quality;
                  if(~iscell(Quality))
                    Quality = {Quality};
                  end
                  
                  n_masks = numel(Quality{1}(1).q);
                  img_2_whole_ids{i} = (previous_id+1):previous_id+n_masks;
                  if(~isempty(img_2_whole_ids{i}))
                    previous_id = img_2_whole_ids{i}(end);
                  end
                  
                  whole_2_img_ids{i} = i*ones(numel(img_2_whole_ids{i}),1);
                  global_2_local_whole_ids{i} = (1:n_masks)';
                  img_overlaps{i} = Quality{1};
                end                                
                
                %bndboxes = obj.get_bboxes_img_ids(1:numel(obj.img_names));
                obj.img_name_to_id_hash = img_name_to_id_hash;

                joint_qualities_file = [obj.exp_dir 'SegmentEval/' obj.mask_type '/' obj.overlap_type '/' obj.imgset '.mat'];
                if(~exist(joint_qualities_file, 'file'))
                  SvmSegm_segment_quality_JL_join(obj.exp_dir, {obj.imgset}, obj.mask_type, obj.overlap_type);
                end
                
                var = load(joint_qualities_file);                
                counter = 1;
                global_2_local_obj_ids = zeros(numel(var.Quality),1);
                for i=1:numel(var.Quality)
                  obj.obj_2_img_ids(i) = obj.img_names_to_ids(var.Quality(i).image_name);
                  if(i~=1 && (obj.obj_2_img_ids(i) ~= obj.obj_2_img_ids(i-1)))
                    counter = 1;
                  end
                  global_2_local_obj_ids(i) = counter;
                  counter = counter + 1;
                end
                
                global_2_local_whole_ids = cell2mat(global_2_local_whole_ids');                 
                obj.whole_2_img_ids = cell2mat(whole_2_img_ids);
                obj.img_2_whole_ids = img_2_whole_ids;
                obj.global_2_local_whole_ids = global_2_local_whole_ids;
                obj.global_2_local_obj_ids = global_2_local_obj_ids;
                obj.img_overlaps = img_overlaps;
                obj_2_img_ids = obj.obj_2_img_ids;
                save(obj.all_whole_ids_file, 'img_name_to_id_hash', 'global_2_local_obj_ids', 'obj_2_img_ids', 'img_2_whole_ids', 'whole_2_img_ids', 'global_2_local_whole_ids', 'img_overlaps', '-V6');
            else
                var = load(obj.all_whole_ids_file);
                obj.img_2_whole_ids = var.img_2_whole_ids;
                obj.whole_2_img_ids = var.whole_2_img_ids;
                if(iscell(obj.whole_2_img_ids))
                    obj.whole_2_img_ids = cell2mat(var.whole_2_img_ids);
                end
                obj.global_2_local_whole_ids = var.global_2_local_whole_ids;
                obj.global_2_local_obj_ids = var.global_2_local_obj_ids;
                obj.img_overlaps = var.img_overlaps;  
                obj.img_name_to_id_hash = var.img_name_to_id_hash;
                obj.obj_2_img_ids = var.obj_2_img_ids;
            end
         end

        function obj = change_exp_dir(obj, new_exp_dir)
          obj.exp_dir = new_exp_dir;          
        end                

        function [Feats, dims] = get_whole_sparse_feats(obj, whole_ids, feat_types, scaling_type)         
          % organize wholes into images
          img_ids = obj.whole_2_img_ids(whole_ids);
          [un_img_ids, ids_a, ids_b] = unique(img_ids);
          
          % load one img to get dimensions
          d = 0;
          for i=1:numel(feat_types)
            D = myload([obj.exp_dir 'MyMeasurements/' obj.mask_type '_' feat_types{i} '/' obj.img_names{1} '.mat'], 'D');
            d = d + size(D,1);
            dims(i) = size(D,1);
          end                                                            
          
          % speed is important here
          mask_type = obj.mask_type;
          img_names = obj.img_names;
          global_2_local_whole_ids = obj.global_2_local_whole_ids;
          exp_dir = obj.exp_dir;

          whole_ranges = [0; find(diff(img_ids)); numel(whole_ids)];
          feat_ranges = [0 cumsum(dims)];
          Feats = [];
          assert(numel(feat_types)==1);
          
          for j=1:numel(feat_types)
            counter = 1;
            range{j} = (feat_ranges(j)+1):feat_ranges(j+1);
            
            ids_x = cell(numel(un_img_ids),1);
            ids_y = cell(numel(un_img_ids),1);
            val = cell(numel(un_img_ids),1);
            n_cols = 0;
            for i=1:numel(un_img_ids)
              these_wholes = whole_ranges(i)+1:whole_ranges(i+1);
              local_ids = global_2_local_whole_ids(whole_ids(these_wholes));
                      
              D = myload([exp_dir 'MyMeasurements/' mask_type '_' feat_types{j} '/' img_names{un_img_ids(i)} '.mat'], 'D');               
              %assert(size(D,2) == numel(obj.img_2_whole_ids{un_img_ids(i)}));    
              
              [ids_x{i},ids_y{i},val{i}] = find(D(:,local_ids));
              ids_y{i} = ids_y{i} + n_cols;
              
              n_cols = n_cols + numel(local_ids);
            end                                 
          end
          
          Feats = sparse(cell2mat(ids_x), cell2mat(ids_y), cell2mat(val), size(D,1), n_cols);
          assert(n_cols == numel(whole_ids));
        end
        
        function [Feats, dims, scaling] = get_whole_feats(obj, whole_ids, feat_types, scaling_type, weights, duh, lut, ndims)
          DefaultVal('*weights', '[]');
          DefaultVal('*lut', 'false');
          DefaultVal('*ndims', '[]');
          
          sorted = sort(whole_ids, 'ascend');
          assert(all(sorted==whole_ids));
          
          % organize wholes into images
          img_ids = obj.whole_2_img_ids(whole_ids);
          [un_img_ids, ids_a, ids_b] = unique(img_ids);
          
          % load one img to get dimensions
          d = 0;

          %the_folder = 'MyMeasurements_ECCV2012/';
          the_folder = 'MyMeasurements/';
            
          for i=1:numel(feat_types)              
            D = myload([obj.exp_dir the_folder obj.mask_type '_' feat_types{i} '/' obj.img_names{1} '.mat'], 'D');
            d = d + size(D,1);
            dims(i) = size(D,1);
          end
                  
          % get features for each image, keep those for the requested wholes         
          Feats = zeros(d, numel(whole_ids), 'single');                                               
          
          % speed is important here
          mask_type = obj.mask_type;
          img_names = obj.img_names;
          global_2_local_whole_ids = obj.global_2_local_whole_ids;          
          exp_dir = obj.exp_dir;

          whole_ranges = [0; find(diff(img_ids)); numel(whole_ids)];
          feat_ranges = [0 cumsum(dims)];
          
          for j=1:numel(feat_types)
            counter = 1;
            range{j} = (feat_ranges(j)+1):feat_ranges(j+1);
            
            for i=1:numel(un_img_ids)
              %vgg_progressbar('feature loading', i/numel(un_img_ids), 5);
                
              these_wholes = whole_ranges(i)+1:whole_ranges(i+1);
              local_ids = global_2_local_whole_ids(whole_ids(these_wholes));
              
              if(~lut)
                D = myload([exp_dir the_folder mask_type '_' feat_types{j} '/' img_names{un_img_ids(i)} '.mat'], 'D');
              else
                D = myload([exp_dir the_folder mask_type '_' feat_types{j} '/' img_names{un_img_ids(i)} '.mat'], 'D');
                lut_max = myload([exp_dir the_folder mask_type '_' feat_types{j} '/' img_names{un_img_ids(i)} '.mat'], 'lut_max');
                lut_min = myload([exp_dir the_folder mask_type '_' feat_types{j} '/' img_names{un_img_ids(i)} '.mat'], 'lut_min');
                
                if(~isempty(ndims))
                  D = D(1:ndims,:);
                end
                
                D = single(D);
                lut_min = lut_min; lut_max = lut_max;
                precs = (lut_max-lut_min)/255;
                precs = 0:precs:1;
                D = precs(D+1);                
              end            
              
              if(~isempty(local_ids))  
                Feats(range{j}, counter:counter+numel(local_ids)-1) = D(:,local_ids);
                counter = counter + numel(local_ids) ;
              end
            end
            
            if(~isempty(scaling_type))
                if(strcmp(scaling_type, 'norm_2') || strcmp(scaling_type, 'norm_1'))
                    chunks = chunkify(1:size(Feats,2), 10);
                    for k=1:numel(chunks) % necessary if data fills up most of memory
                        [Feats(range{j}, chunks{k}), scaling{j}] = scale_data(Feats(range{j},chunks{k}), scaling_type);
                    end
                else
                    if(~strcmp(scaling_type, 'none'))
                        error('not ready for this');
                    end
                end
            else
              scaling{j} = [];
            end
                        
            if(~isempty(weights) && weights(j)~=1)
                Feats(range{j},:) = Feats(range{j},:)*weights(j);
            end
          end
          
          if((numel(feat_types)>1) && (strcmp(scaling_type, 'norm_2') || strcmp(scaling_type, 'norm_1')))
            chunks = chunkify(1:size(Feats,2), 10);
            for k=1:numel(chunks) % necessary if data fills up most of memory
                Feats(:,chunks{k}) = scale_data(Feats(:,chunks{k}), scaling_type);
            end
          end          
        end
        
        function [overlap] = get_overlap_wholes_objclass(obj, class_id)
          [img_ids, obj_ids, q] = obj.collect_category_imgs(class_id);

          overlap = zeros(numel(obj.whole_2_img_ids),1);
          for i=1:numel(img_ids)
            wholes = cell2mat(obj.img_2_whole_ids(img_ids(i)));

            for j=1:numel(q{i})
              overlap(wholes) = max(overlap(wholes), q{i}(j).q);
            end
          end
        end
        
        function obj_ids = whole_2_obj_ids(obj, whole_ids)
            img_ids = obj.whole_2_img_ids(whole_ids);
            obj_ids = zeros(numel(whole_ids),1);
            img_overlaps = obj.img_overlaps;
            obj_2_img_ids = obj.obj_2_img_ids;
            global_2_local_whole_ids = obj.global_2_local_whole_ids;
            parfor j=1:numel(whole_ids)                
                img_ov = img_overlaps{img_ids(j)};
                global_obj_ids = find(obj_2_img_ids==img_ids(j));
                q = [img_ov(:).q];
                [q, pos] = max(q(global_2_local_whole_ids(whole_ids(j)), :));
                if(q>0)
                    obj_ids(j) = global_obj_ids(pos);
                else
                    obj_ids(j) = -1;
                end
            end
        end
        
        function [overlaps, labels] = get_overlaps_img(obj, img_name)
          if(~isstr(img_name))
              img_name = obj.img_names{img_name};
          end
          
          load([obj.exp_dir 'SegmentEval/' obj.mask_type '/' obj.overlap_type '/' img_name '.mat']);
          overlaps = [Quality{1}.q];
          labels = [Quality{1}.class];
          %[val, id] = max(overlap)
          %real_id = obj.show_best_masks(img_name);
        end

        function [overlap] = get_overlaps_wholes(obj, whole_ids)    
          overlap = zeros(numel(obj.whole_2_img_ids),numel(obj.categories));
          img_ids = obj.whole_2_img_ids;
          un_img_ids = unique(img_ids);
          
          img_overlaps = obj.img_overlaps;
          counter = 1;
          for i=1:numel(un_img_ids)
            %var = load([obj.exp_dir 'SegmentEval/' obj.mask_type '/' obj.overlap_type '/' obj.img_names{un_img_ids(i)} '.mat']);
            %Q = var.Quality;            
            %q = [Q{1}.q];
            %classes = [Q{1}.class];            
            duh = img_overlaps{i};
            q = [duh.q];
            classes = [duh.class];
            [un_classes, ids] = unique(classes);
            n_segms = size(q,1);
            n_un_classes = numel(un_classes);
            
            newq = zeros(n_segms, n_un_classes);
            
            for j=1:n_un_classes
              newq(:,j) = max(q(:, classes == un_classes(j)), [], 2);
            end
            
            counter2 = counter + size(q,1) - 1;
            overlap(counter:counter2, un_classes) = newq;
            
            counter = counter2 + 1;
          end
          
          overlap = overlap(whole_ids,:);
        end        
                
        function [trunc, occl] = are_wholes_trunc(obj, whole_ids)
          % returns whether the objects having largest bbox overlap with segment bbox
          % are marked as truncated on the detection annotations
          trunc = false(numel(whole_ids),1);
          occl = trunc;
          [info] = obj.get_info_bbox_det(obj.whole_2_img_ids(whole_ids));
          local_whole_ids = obj.global_2_local_whole_ids(whole_ids);
          for i=1:numel(whole_ids)
            Q = [info{i}(:).q];
            [~, obj_id] = max(Q,[],2);
            obj_id = obj_id(local_whole_ids(i));
            trunc(i) = info{i}(obj_id).truncated;
            occl(i) = info{i}(obj_id).occluded;
          end
        end
        
        function [info] = get_info_bbox_det(obj, img_ids)
          %%%% retrieves overlaps between segment bounding boxes and
          %%%% detection challenge annotated bounding boxes as well as
          %%%% annotations (e.g. truncation)
          cache = [obj.exp_dir 'SegmentEval/' obj.mask_type '/' obj.imgset '_bndboxes_det_gt.mat'];
          
          if(exist(cache,'file'))
            var = load(cache);
            info = var.Quality;
            info = info(img_ids);
            return;
          else      
            % get segment bounding boxes
            bndboxes = obj.get_bboxes_img_ids(1:numel(obj.img_names));
            
            npos=0;
            gtids = obj.img_names;
            gt(length(gtids))=struct('BB',[],'diff',[],'det',[]);
            Quality = cell(numel(gtids),1);
            for i=1:length(gtids)
              try 
                  rec=PASreadrecord(sprintf([obj.exp_dir 'Annotations/%s.xml'],gtids{i}));
              catch
                  % (dummy for our synthetic single-object-per-image dataset)
                  rec = [];
                  rec.objects(1).class = 'aeroplane';
                  rec.objects(1).truncated = false;
                  rec.objects(1).occluded = false;
                  rec.objects(1).view = [];
                  rec.objects(1).label = -1;
                  rec.objects(1).bbox = [1 1 2 2];
                  rec.objects(1).difficult = false;
              end            
              
              bbgt=cat(1,rec.objects.bbox)';
              diff=[rec.objects.difficult];
              
              for j=1:numel(rec.objects)                
                %bbgt = bbgt';
                bb = bndboxes{i};
                bi=[max(bb(:,1),bbgt(1,j)) max(bb(:,2),bbgt(2,j)) min(bb(:,3),bbgt(3,j)) min(bb(:,4),bbgt(4,j))];
                q = zeros(size(bndboxes{i},1),1);
                iw=bi(:,3)-bi(:,1)+1;
                ih=bi(:,4)-bi(:,2)+1;
                nz = iw>0 & ih>0;
                
                % compute overlap as area of intersection / area of union
                ua=(bb(:,3)-bb(:,1)+1).*(bb(:,4)-bb(:,2)+1)+...
                  (bbgt(3,j)-bbgt(1,j)+1)*(bbgt(4,j)-bbgt(2,j)+1)-...
                  iw.*ih;
                ov=(iw.*ih)./ua;
                q(nz) = ov(nz);
                
                Quality{i}(j).class = VOC09_classname_to_id(rec.objects(j).class);
                Quality{i}(j).object = j;
                Quality{i}(j).q = single(q);
                Quality{i}(j).diff = diff(j);
                Quality{i}(j).truncated = rec.objects(j).truncated;
                Quality{i}(j).occluded = rec.objects(j).occluded;
                Quality{i}(j).view = rec.objects(j).view;
                Quality{i}(j).label = rec.objects(j).label;


              end
            end
            save(cache, 'Quality', '-V6');
            info = Quality;
            info= info(img_ids);
          end
        end
        
        function [labels] = get_labels_wholes(obj)
            overlaps = zeros(numel(obj.whole_2_img_ids),numel(obj.categories));
            for i=1:numel(obj.categories)
                [overlaps(:,i)] = get_overlap_wholes_objclass(obj, i);
            end
            [val, labels] = max(overlaps');
        end

        function out = whole_fun(obj, fun)
            img_names = obj.img_names(obj.whole_2_img_ids);
            local_whole_ids = obj.global_2_local_whole_ids(1:numel(img_names));

            [un_imgs, n] = unique(img_names); 
            n = [0; n];
            counter = 1;     
            out = cell(numel(img_names),1);
            mask_type = obj.mask_type;
            exp_dir = obj.exp_dir;
            
            t = tic();
            parfor i=1:numel(un_imgs)
            %for i=1:numel(un_imgs)
                i
                masks_wholes = myload([exp_dir 'MySegmentsMat/' mask_type '/' un_imgs{i} '.mat'], 'masks');
                                                
                if(iscell(masks_wholes))
                    masks_wholes = masks_wholes{1};
                end

                img_mask_base_id = n(i)+1;
                for  j=1:(n(i+1)-img_mask_base_id+1)
                    out{i}(j,:) = fun(masks_wholes(:,:,j));
                end       
                %n(1) = [];
            end           
            toc(t)
        end

        function aspect_ratio = get_aspect_ratio_wholes(obj)
            function aratio = asp_ratio(mask)
                [i,j] = find(mask);
                i_len = max(i) - min(i);
                j_len = max(j) - min(j);
                aratio = i_len/j_len;
            end

            aspect_ratio = cell2mat(obj.whole_fun(@asp_ratio));
        end 
        
        function bbox = get_norm_bbox_wholes(obj)
            function bbox = get_norm_bbox(mask)
                [i,j] = find(mask);
                bbox = [min(i) min(j) max(i) max(j)];
                bbox([1 3]) = bbox([1 3])/size(mask,1);
                bbox([2 4]) = bbox([2 4])/size(mask,2);
            end

            bbox = cell2mat(obj.whole_fun(@get_norm_bbox));
        end                

        function [best_wholes_ids, q] = get_best_wholes(obj, obj_ids)         
            var = load([obj.exp_dir 'SegmentEval/' obj.mask_type '/' obj.overlap_type '/' obj.imgset '.mat']);
            Quality = var.Quality;
            best_wholes_ids = zeros(numel(obj_ids),1);
            q = best_wholes_ids;
            for i=1:numel(obj_ids)
                [q(i), best_local_whole_id] = max(Quality(obj_ids(i)).q);
                whole_ids = obj.img_2_whole_ids(obj.obj_2_img_ids(obj_ids(i)));
                best_wholes_ids(i) = whole_ids{1}(best_local_whole_id);
            end
            
            [best_wholes_ids, srt] = sort(best_wholes_ids, 'ascend');
            q = q(srt);
        end
        
        function q = get_img_mask_quality(obj, img_name)
          var = load([obj.exp_dir 'SegmentEval/' obj.mask_type '/' obj.overlap_type '/' img_name '.mat']);
          q = var.Quality;
          if(numel(q) == 1)
            q = q{1};
          end
        end

        function [best_whole_ids, best_whole_q] = get_best_mask_ids(obj, img_name)    
           best_whole_ids = [];
           best_whole_q = [];
           if(~iscell(img_name))
            img_name = obj.img_names{img_name};
           end             
           var = load([obj.exp_dir 'SegmentEval/' obj.mask_type '/' obj.overlap_type '/' img_name '.mat']);
           Q = var.Quality{1};
           % Q = obj.img_overlaps{img_id}; % this is really slow, for some
           % reason
           
           for i=1:numel(Q)
            [best_whole_q(i), best_whole_ids(i)] = max(Q(i).q);
           end
        end           
        
        function [whole_ids] = get_whole_ids_in_qual_range(obj, categ, max_qual, min_qual)
          DefaultVal('*min_qual', '-inf');
          DefaultVal('*max_qual', 'inf');
          DefaultVal('*categ', '[]');
                    
          whole_ids = [];
          q_imgs = cellfun(@(a) [a.q], obj.img_overlaps, 'UniformOutput', false);
          larger = cell(numel(q_imgs),1);
          smaller = larger;
          for i=1:numel(q_imgs)
            max_q = max(q_imgs{i}, [], 2);
            larger{i} = max_q>min_qual;
            smaller{i} = max_qual> max_q;
            
            passing = larger{i} & smaller{i};
            whole_ids = [whole_ids; obj.img_2_whole_ids{i}(passing)'];
          end                     
        end
          
        function [img_ids, obj_ids, Q] = collect_category_imgs(obj, category)   
            if(~isstr(category))
                category = obj.categories{category};
            end

            var = load([obj.exp_dir 'SegmentEval/' obj.mask_type '/' obj.overlap_type '/' obj.imgset '.mat']);
            all_Quality = var.Quality;

            imgs = {};
            obj_ids = [];

            for i=1:numel(all_Quality)
                if(all_Quality(i).class == VOC09_classname_to_id(category))
                    imgs = [imgs all_Quality(i).image_name];
                    obj_ids = [obj_ids; i];
                end
            end
            
            [imgs, a, b] = unique(imgs);

            Q = cell(numel(imgs),1);
            new_obj_ids = Q;
            for i=1:numel(imgs)                           
                new_obj_ids{i} = [obj_ids(b==i)];
                Q{i} = all_Quality(new_obj_ids{i});
            end

            obj_ids = new_obj_ids;
            img_ids = obj.img_names_to_ids(imgs);
        end

        function best_whole_ids = get_best_whole_ids(obj, obj_ids)
            assert(~iscell(obj_ids));
            var = load([obj.exp_dir 'SegmentEval/' obj.mask_type '/' obj.overlap_type '/' obj.imgset '.mat']);
            all_Quality = var.Quality;

            first_segm_id_obj_id = zeros(numel(all_Quality),1);
            first_segm_id_obj_id(1) = 1;
            img_name = all_Quality(1).image_name;
            for i=2:numel(all_Quality)
                %all_Quality(i).image_name
                if(~strcmp(all_Quality(i).image_name, img_name))
                    first_segm_id_obj_id(i) = first_segm_id_obj_id(i-1) + numel(all_Quality(i-1).q);
                else
                    first_segm_id_obj_id(i) = first_segm_id_obj_id(i-1);
                end
                
                img_name = all_Quality(i).image_name;
            end
            
            best_whole_ids = zeros(numel(obj_ids),1);
            for i=1:numel(obj_ids)
                [val, id] = max(all_Quality(obj_ids(i)).q);
                best_whole_ids(i) = first_segm_id_obj_id(obj_ids(i)) + id -1;
            end
        end        

        function img_ids = img_names_to_ids(obj, img_names)
            if(~iscell(img_names))
                img_names = {img_names};
            end
            
            img_ids = zeros(numel(img_names),1);
            for j=1:numel(img_names)
              %t1 = tic();
              
              %%%% this was here %%%%
              %val = str2double(img_names{j}([1:4 6:end]));
              %img_ids(j) = find(ismember(obj.img_name_to_id_hash, val));
              %%%%%%%%%%%%%%%%%%%%%%%
              % toc(t1)
              
              img_ids(j) = find(strcmp(obj.img_names, img_names{j}));
              %toc(t2)
            end
        end

        function part_ids = get_part_ids(obj, segm_ids)            
            for i=1:numel(segm_ids)
                part_ids{i} = obj.whole_2_part_ids{segm_ids(i)};
            end            
        end

        function [img_names, local_whole_ids, local_part_ids] = get_part_id_origins(obj,part_ids)
            whole_ids = obj.part_2_whole_ids(part_ids);            
            local_whole_ids = obj.global_2_local_whole_ids(whole_ids);            
            local_part_ids = obj.global_2_local_part_ids(part_ids);
            img_names = obj.img_names(obj.whole_2_img_ids(whole_ids));
        end

        function [Imgs, masks_wholes_to_show] = get_whole_Imgs_and_masks(obj, whole_ids)
            Imgs = obj.get_whole_Imgs(whole_ids);
            masks_wholes_to_show = obj.get_masks(whole_ids, 'whole');
        end

        function [Imgs, tl_bb_corners] = get_object_Imgs(obj, obj_ids, just_bbox)
            DefaultVal('*just_bbox', 'false');
            img_names = obj.img_names(obj.obj_2_img_ids(obj_ids));            
            Imgs = obj.get_Imgs(img_names);
            if(just_bbox)
                [Imgs, tl_bb_corners] = obj.cut_bbox(Imgs, obj_ids);
            else
                tl_bb_corners = zeros(2, numel(obj_ids));
            end
        end

        function [masks, sp_app, sp] = get_img_masks(obj, img_id)
            var = load([obj.exp_dir 'MySegmentsMat/' obj.mask_type '/' obj.img_names{img_id} '.mat']);
            
            masks = var.masks;
            if(isfield(var, 'sp_app'))
                % superpixel approximation
                sp_app = var.sp_app;
                % superpixel map
                sp = var.sp;
            else
                sp_app = [];
                sp = [];
            end
        end
        
        function Imgs = get_whole_Imgs(obj, whole_ids)
            img_names = obj.img_names(obj.whole_2_img_ids(whole_ids));            
            Imgs = get_Imgs(obj, img_names);
        end

        function Imgs = get_Imgs(obj, img_names)
          if(~iscell(img_names))
            img_names = obj.img_names(img_names);
          end
          
          % if they don't come sorted it doesn't work!
          %[un_imgs, n] = unique(img_names);
          un_imgs = img_names;
          n = 1:numel(un_imgs);
          
          Imgs = cell(numel(img_names),1);
          counter = 1;
          for i=1:numel(un_imgs);
            nextI = imread([obj.exp_dir 'JPEGImages/' un_imgs{i} '.jpg']);
            for  j=counter:n(1)
              Imgs{counter} = nextI;
              counter = counter + 1;
            end
            n(1) = [];
          end
        end
        
        function [bbox_Imgs, tl_bb_corners] = cut_bbox(obj, Imgs, obj_ids)            
            bbox_Imgs = cell(numel(Imgs),1);
            tl_bb_corners = zeros(2, numel(obj_ids));
            
            img_ids = unique(obj.obj_2_img_ids);
            
            if(strcmp(obj.overlap_type, 'overlap')) % segmentation gt
                bndboxes = obj.get_bboxes_img_ids(img_ids, 'ground_truth');                    
                gt_bndboxes = cell2mat(bndboxes);
                bndbox = gt_bndboxes(obj_ids,:);
            elseif(strcmp(obj.overlap_type, 'overlap_and_bbox_edges')) % bbox gt
                bndbox = obj.get_bboxes(obj_ids);    
                bndbox = [[bndbox(:).xmin]' [bndbox(:).ymin]' [bndbox(:).xmax]' [bndbox(:).ymax]'];
            end                            
            
            %bndbox = obj.get_bboxes(obj_ids);            
            for i=1:numel(obj_ids)
                %tl_bb_corners(:,i) = [bndbox(i).xmin bndbox(i).ymin];
                %range_x = bndbox(i).xmin:bndbox(i).xmax;
                %range_y = bndbox(i).ymin:bndbox(i).ymax;
                
                tl_bb_corners(:,i) = bndbox(i, [1 2]);
                range_x = bndbox(i,1):bndbox(i,3)-1;
                range_y = bndbox(i,2):bndbox(i,4)-1;
                bbox_Imgs{i} = Imgs{i}(range_y,range_x,:);
            end
        end
        
        function bndboxes = get_bboxes(obj, obj_ids)                    
            img_names = obj.img_names(obj.obj_2_img_ids(obj_ids));
            for i=1:numel(obj_ids)
                rec=PASreadrecord(sprintf(obj.VOCopts.annopath,img_names{i}));
                bndboxes(i) = rec.objects(obj.global_2_local_obj_ids(obj_ids(i))).bndbox;
            end            
        end
        
        function bndboxes = get_bboxes_from_masks(obj, masks)
            % assumes one mask per cell
            if(iscell(masks))
                for i=1:numel(masks)
                    [y,x] = find(masks{i});
                    bndboxes(i).xmax = max(x);
                    bndboxes(i).xmin = min(x);
                    bndboxes(i).ymax = max(y);
                    bndboxes(i).ymin = min(y);
                end
            else 
                error('to implement later');
            end
        end
        
        function bndboxes = get_bboxes_img_ids(obj, img_ids, other_mask_type)
          % pascal format bboxes (x1 y1 x2 y2) from all segments
          bndboxes = cell(numel(img_ids),1);          
          
          if(nargin == 2)
            cache = [obj.exp_dir 'MySegmentsMat/' obj.mask_type '/' obj.imgset '_bndboxes.mat'];
            browser = obj;
          else
            cache = [obj.exp_dir 'MySegmentsMat/' other_mask_type '/' obj.imgset '_bndboxes.mat'];
            browser = SegmBrowser(obj.exp_dir, other_mask_type, obj.imgset);
          end
          
          if(exist(cache, 'file'))
            load(cache);
          else
            disp('Caching bounding box overlaps.');
            parfor i=1:numel(obj.img_names)
            %for i=1:numel(obj.img_names)            
              masks = browser.get_img_masks(i);
              if(iscell(masks))
                  masks = masks{1};
              end
              
              bndboxes{i} = zeros(size(masks,3), 4);
              for j=1:size(masks,3)
                bndboxes{i}(j,:) = regionprops_BB_mine(masks(:,:,j));
                bndboxes{i}(j,3) = bndboxes{i}(j,1) + bndboxes{i}(j,3);
                bndboxes{i}(j,4) = bndboxes{i}(j,2) + bndboxes{i}(j,4);
              end
            end
            save(cache, 'bndboxes', '-V6');
          end
          
          bndboxes = bndboxes(img_ids);
        end          
        
        function masks_to_show = get_masks(obj, ids, whole_or_part)
            DefaultVal('*whole_or_part', '''whole''');
            
            is_part = true;
            if(strcmp(whole_or_part, 'whole'))
                is_part = false;
            end
            
            if(~is_part)
                img_names = obj.img_names(obj.whole_2_img_ids(ids));                
                local_ids = obj.global_2_local_whole_ids(ids);
                inner_ids = ones(numel(local_ids),1);
            elseif(is_part)
                [~, inner_ids, local_ids] = obj.get_part_id_origins(ids);
                %local_ids = obj.global_2_local_part_ids(ids);
                img_names = obj.img_names(obj.whole_2_img_ids(obj.part_2_whole_ids(ids)));
            else
                error('no such type');
            end            
                        
            masks_to_show = cell(numel(img_names),1);
            un_imgs = img_names; n = 1:numel(img_names);
            %[un_imgs, n] = unique(img_names);
            counter = 1;
            for i=1:numel(un_imgs)
                if(is_part)
                    var = load([obj.exp_dir 'MyPartSegmentsMat/' obj.mask_type '/' un_imgs{i} '.mat'], 'masks');
                else
                    var = load([obj.exp_dir 'MySegmentsMat/' obj.mask_type '/' un_imgs{i} '.mat'], 'masks');
                end                
                                
                masks = var.masks;                
                if(~iscell(masks))
                  masks = {masks};
                end
                
                for  k=counter:n(1)
                    masks_to_show{counter} = masks{inner_ids(counter)}(:,:,local_ids(counter));
                    counter = counter + 1;
                end
                
                n(1) = [];
                
                if 0 %%%% debug %%%%
                    n_masks = size(masks,3)
                    n(1)-counter+1
                    subplot(1,2,1), obj.show_imgs(un_imgs(i));
                end
            end

            %masks_to_show = masks_to_show(order);
        end
          
        function reg_masks = get_registered_masks(obj, ids, max_height, max_width)          
          DefaultVal('*max_width', '[]');
          img_names = obj.img_names(obj.whole_2_img_ids(ids));
          local_ids = obj.global_2_local_whole_ids(ids);
          inner_ids = ones(numel(local_ids),1);
       
          un_imgs = img_names; n = 1:numel(img_names);
          %[un_imgs, n] = unique(img_names);
          counter = 1;
          for i=1:numel(un_imgs)
            i
            if(i==1 || ~strcmp(un_imgs{i}, previous_img))
              var = load([obj.exp_dir 'MySegmentsMat/' obj.mask_type '/' un_imgs{i} '.mat'], 'masks');            
              masks = var.masks;
              if(~iscell(masks))
                masks = {masks};
              end
              previous_img = un_imgs{i};
            end            
            for  k=counter:n(1)
              mask = masks{inner_ids(counter)}(:,:,local_ids(counter));
              bbox = regionprops_BB_mine(mask);
              mask2 = mask(bbox(2):bbox(2)+bbox(4)-1, bbox(1):bbox(1)+bbox(3)-1);
              height_ratio = size(mask2, 1)/max_height;
              des_size = floor(size(mask2)/height_ratio);              
              des_size(1) = max_height;
              
              mask3 = imresize(mask2,des_size);
                            
              %reg_masks{counter} = mask3;
              pre_reg_masks{counter} = logical(mask3);
              counter = counter + 1;
            end
                        
            n(1) = [];
          end
                    
          if(isempty(max_width))
            max_width = max(cellfun(@(a) size(a,2), pre_reg_masks));                  
          end
          
          reg_masks = false(max_height, max_width, numel(ids));
          for i=1:numel(un_imgs)
            extra = mod(max_width - size(pre_reg_masks{i},2), 2);
            div = floor((max_width - size(pre_reg_masks{i},2))/2);
            
            padsize_pre = [0 div];
            tmp_mask = padarray(pre_reg_masks{i}, padsize_pre, false, 'pre');
            padsize_post = [0 div+extra];
            
            tmp_mask = padarray(tmp_mask, padsize_post, false, 'post');
            reg_masks(:,:,i) = tmp_mask;
          end
        end
        
        %%%% Visualization %%%%

        function Imgs = show_objects(obj, obj_ids)
            just_bbox = true;
            [Imgs, tl_bb_corners] = obj.get_object_Imgs(obj_ids, just_bbox);
            if nargout == 0
                montage_new(Imgs);
            end
        end

        function show_imgs(obj, img_ids)
            false;
            if(iscell(img_ids(1)))
                [Imgs] = obj.get_Imgs(img_ids); % it has names, strings
            else            
                [Imgs] = obj.get_Imgs(obj.img_names(img_ids));
            end
            montage_new(Imgs);
        end

        function Imgs = show_wholes(obj, whole_ids, titles)          
          DefaultVal('*titles', '[]');
          masks_wholes_to_show = obj.get_masks(whole_ids, 'whole')
          Imgs = obj.get_Imgs(obj.img_names(obj.whole_2_img_ids(whole_ids)));
            
          Imgs = subplot_auto_transparent_multiple_imgs(masks_wholes_to_show, Imgs, titles);
        end
               
        function best_id = show_best_masks(obj, img_name, n_best)
          if(~isstr(img_name))
            img_name = obj.img_names{img_name};
          end
          
          DefaultVal('*n_best', '1');
          var = load([obj.exp_dir 'MySegmentsMat/' obj.mask_type '/' img_name '.mat']);
          var2 = load([obj.exp_dir 'SegmentEval/' obj.mask_type '/' obj.overlap_type '/' img_name '.mat']);
          masks = var.masks;
          Q = var2.Quality;
          I = imread([obj.exp_dir 'JPEGImages/' img_name '.jpg']);
          [imgtoshow, score, best_id] = SvmSegm_show_best_segments(I, Q, masks, n_best);
          max(score)
        end
        
        function categ_masks = collect_all_category_best_masks(obj, category, img_names)
          DefaultVal('*img_names', '[]');
          if(isempty(img_names))            
            [img_ids, obj_ids, Q] = obj.collect_category_imgs(category);
            %img_names = obj.img_names;
          end
          
          VOCinit();
          
          img_names = obj.img_names(img_ids);
          
          counter = 1;          
          for i=1:numel(img_ids)                                    
            whole_ids = obj.get_best_wholes(obj_ids{i});
            var = load([obj.exp_dir 'MySegmentsMat/' obj.mask_type '/' img_names{i} '.mat']);
            masks = var.masks;
            
            %var2 = load([obj.exp_dir 'SegmentEval/' obj.mask_type '/' obj.overlap_type '/' img_names{i} '.mat']);            
            %Q = var2.Quality;
            local_ids = obj.global_2_local_whole_ids(whole_ids);
            for j=1:numel(local_ids)
              categ_masks{counter} = masks(:,:,local_ids(j));
              counter = counter + 1;
            end
          end
        end
        
        function show_all_category_imgs(obj, category, img_names)
          DefaultVal('*img_names', '[]');
          if(isempty(img_names))
            img_names = obj.img_names;
          end
          
          VOCinit();
          
          imgs = {};
          counter = 1;
          for i=1:numel(img_names)
            rec(i)=PASreadrecord(sprintf(VOCopts.annopath,img_names{i}));
            
            img_loaded = false;
            for j=1:numel(rec(i).objects);
              if(strcmp(rec(i).objects(j).class, category))
                if(~img_loaded)
                  img = imread([obj.exp_dir 'JPEGImages/' img_names{i} '.jpg']);
                  img_loaded = true;
                end
                % get only region of interest
                range_x = rec(i).objects(j).bndbox.xmin:rec(i).objects(j).bndbox.xmax;
                range_y = rec(i).objects(j).bndbox.ymin:rec(i).objects(j).bndbox.ymax;
                imgs{counter} = img(range_y, range_x,:);
                counter = counter + 1;
              end
            end
          end
          
          montage_new(imgs);
        end                        
        
        function voc_segm_visuals(obj, png_folder, jpeg_folder, name, plot_labels)  
            DefaultVal('*plot_labels', 'true');
            dir = [obj.exp_dir '/visuals/' name '_visuals/'];
            if(~exist(dir))
                mkdir(dir);
            end

            for i=1:numel(obj.img_names)
                Ipng = imread([png_folder obj.img_names{i} '.png']);
                Ijpg = imread([jpeg_folder obj.img_names{i} '.jpg']);
                
                labels = unique(Ipng);
                labels(labels==0) = [];
                
                masks = false(size(Ipng,1), size(Ipng,2), numel(labels));                
                for j=1:numel(labels)
                    masks(:,:,j) = (Ipng==labels(j));
                end
                
                cmap = VOClabelcolormap(256);
              
                
                if 1
                    intens = ones(numel(labels),1);
                    transp = 0.7;
                    newI = create_transparent_multiple_colors(masks, Ijpg, true, labels, transp, intens);
                    voc_labels = labels;

                    line_width = 2;                    
                    edge_I = zeros(size(masks));
                    one_img_edge_I = zeros([size(Ijpg,1) size(Ijpg,2)]);
                    for j=1:size(masks,3)
                        b{j} = cell2mat(bwboundaries(masks(:,:,j)));
                        this_edge_I = zeros([size(Ijpg,1) size(Ijpg,2)]);
                        
                        ind = sub2ind([size(Ijpg,1) size(Ijpg,2)], b{j}(:,1), b{j}(:,2));
                        this_edge_I(ind) = 1;
                        
                        if(line_width==0)
                            this_edge_I = logical(this_edge_I);
                        elseif(line_width==1) % beautiful
                            this_edge_I = bwmorph(this_edge_I, 'dilate');
                        else
                            this_edge_I = bwmorph(this_edge_I, 'dilate');
                            this_edge_I = bwmorph(this_edge_I, 'dilate');
                        end
                        
                        edge_I(:,:,j) = this_edge_I;
                        one_img_edge_I(this_edge_I) = j+1;
                    end                    
                    
                    if(~isempty(voc_labels))
                        for z=1:numel(voc_labels)
                            one_img_edge_I(one_img_edge_I==(z+1)) = voc_labels(z)+1;
                        end
                    end
                    one_img_edge_I= ind2rgb(one_img_edge_I, cmap);
                    one_img_edge_I = uint8(one_img_edge_I*255);
                    
                    for j=1:3
                        newI_j = newI(:,:,j);
                        one_img_j = one_img_edge_I(:,:,j);
                        to_change = sum(one_img_edge_I,3)~=0;
                        newI_j(to_change) = one_img_j(to_change);
                        newI(:,:, j) = newI_j;
                    end                                        
                else                    
                    newI = make_plot_segments_with_overlap(masks, Ijpg, Ipng, labels, cmap);            
                end
                
                sc(newI);
                this_cmap = ones(size(masks,3), 3);
                if(plot_labels)
                    SvmSegm_plot_labels_at_segms(masks, obj.categories(labels), this_cmap)
                end
                duh = getframe();
                imwrite(duh.cdata, [dir obj.img_names{i} '.png']);
            end
        end
        
        % Dataset consistency
        
        function verify_GT_consistency(obj)
          % checks if 
          n_gt = cellfun(@numel, obj.img_overlaps);                        
          for i=1:numel(obj.img_names)
            retval = vgg_progressbar('Verifying ground truth and qualities consistency. ', i/numel(obj.img_names));
            I = imread([obj.exp_dir 'SegmentationObject/' obj.img_names{i} '.png']);
            u = unique(I);
            u(u==255) = []; % don't care pixels
            u(u==0) = []; % background
                        
            % check that there are as many objects in obj labeled images as
            % qualities
            
            %%% check if there are as many ground truth masks as qualities
            var = load([obj.exp_dir 'MySegmentsMat/ground_truth/' obj.img_names{i}]);
            if(numel(u) ~= size(var.masks,3))
              fprintf('Found inconsistent ground truth masks: %s\n', obj.img_names{i});
              obj.show_imgs(obj.img_names_to_ids(obj.img_names{i}));
              r = input('Would you like to fix gt masks in this image ? I''ll also delete any features extracted on them. (y/n)','s');
              if(strcmp(r, 'y'))
                delete([obj.exp_dir 'MySegmentsMat/ground_truth/' obj.img_names{i} '.mat']);
                system(['rm ' obj.exp_dir 'MyMeasurements/ground_truth_*/' obj.img_names{i} '.mat']);
                SvmSegm_generate_gt_masks(obj.exp_dir, obj.img_names(i), false);
              end              
            end
            
            if(numel(u)~=n_gt(i))
              fprintf('Found inconsistent segment qualities: %s\n', obj.img_names{i});              
              obj.show_imgs(obj.img_names_to_ids(obj.img_names{i}));
              r = input('Would you like to fix segment qualities in this image ? (y/n)','s');
              if(strcmp(r, 'y'))
                delete([obj.exp_dir 'SegmentEval/' obj.mask_type '/' obj.img_names{i} '.mat']);
                SvmSegm_segment_quality_all_JL(obj.img_names(i), obj.exp_dir, obj.mask_type, obj.overlap_type);
              else
                fprintf('\nSkipping!\n');
              end
            end
          end
          
          disp('Ground truth SegmentationObj files and segment qualities are consistent!');
        end        
        
        function problem = detect_corrupted_feat_files(obj, feat_type)
          img_names = obj.img_names;
          exp_dir = obj.exp_dir;
          mask_type = obj.mask_type;
          
 				  problem = {};
          for i=1:numel(img_names)
            try
              load([exp_dir 'MyMeasurements/' mask_type '_' feat_type '/' img_names{i} '.mat'], 'D');
            catch
              disp('problem!');
              problem = cat(2, problem, img_names{i})
							keyboard;
            end
          end
        end
        
        % Features
        
        function all_D = sample_really_dense_sift_imgs(obj, n_per_img, max_N)
          n_imgs = max_N/n_per_img;
          assert(numel(obj.img_names) >= n_imgs);
          
          r = randperm(numel(obj.img_names));
          
          sel_imgs = obj.img_names(r(1:n_imgs));
          all_D = cell(1,numel(sel_imgs));
          parfor i=1:numel(sel_imgs)
            t = tic();
            %vgg_progressbar('sampling sifts', i/numel(sel_imgs), 0.5)
            I = obj.get_Imgs(obj.img_names_to_ids(sel_imgs{i}));
            [F,D] = vl_phow(im2single(I{1}), 'Color', 'gray');
            r = randperm(size(D,2));
            all_D{i} = D(:,r(1:n_per_img));
            toc(t)
          end
          
          all_D = cell2mat(all_D);
        end               

        % Evaluation

        function [ap,class_fp, class_tp] = voc_detection_score(obj, local_ids, global_ids, labels, scores)          
          boxes = obj.get_bboxes_img_ids(1:numel(obj.img_names));          
          assert(numel(local_ids) == numel(boxes));

          class_fp = cell(20,1);
          class_tp = cell(20,1);

          bbox_overlaps = obj.get_info_bbox_det(1:numel(obj.img_names));
          %for i=1:numel(bbox_overlaps)
          %  for j=1:numel(bbox_overlaps{i})
          %    bbox_overlaps{i}(j).q = bbox_overlaps{i}(j).q(ids{i});
          %  end
          %end

          for i=1:20                  
          %parfor i=1:20    
            class_bbox_overlaps = bbox_overlaps;
            %boxes1 = cell(numel(boxes),1);

            boxes1 = cell(numel(boxes),1);
            in_class_global_ids = boxes1;
            for j=1:numel(boxes)    
              labels_right = labels{j}==i;
              boxes1{j} = [boxes{j}(local_ids{j},:) scores{j}'];
              boxes1{j} = boxes1{j}(labels_right,:);           
              in_class_global_ids{j} = global_ids{j}(labels_right);
            end

            for j=1:numel(boxes)              
              to_remove = false(numel(class_bbox_overlaps{j}),1);
              labels_right = labels{j}==i;
              the_I = local_ids{j}(labels_right);
              for k=1:numel(class_bbox_overlaps{j})
                class_bbox_overlaps{j}(k).q = class_bbox_overlaps{j}(k).q(the_I);
                if ((class_bbox_overlaps{j}(k).class~=i))
                  to_remove(k) = true;
                %elseif(isempty(class_bbox_overlaps{j}(k).q))
                  %to_remove(k) = true;
                end
              end

              class_bbox_overlaps{j}(to_remove) = [];
            end                  

            %%% similar to original voc function %%%%
            %[ap(i), fp, tp] = voc_ap(obj.categories{i}, boxes1, obj.imgset);                                    
            %%% fast version %%%
            [ap(i), fp, tp] = voc_ap_fast_and_stupid(boxes1, class_bbox_overlaps); 

            empty = cellfun(@isempty, in_class_global_ids);            
            in_class_global_ids(empty) = [];
            in_class_global_ids = cell2mat(in_class_global_ids');
            class_fp{i} = in_class_global_ids(logical(fp));
            class_tp{i} = in_class_global_ids(logical(tp));
          end          
        end

        function [ap,class_fp, class_tp] = voc_cls_score(obj, local_ids, global_ids, labels, scores)                     
          class_fp = cell(20,1);
          class_tp = cell(20,1);                               
          
          for i=1:20
            class_scores = scores;
            id_winner = zeros(numel(scores),1);
            for j=1:numel(scores)
              % get the max score
              if(~isempty(scores{j}))
                ids = find(labels{j}==i);
                if(any(ids))
                  [class_scores{j}, id_winner(j)] = max(scores{j}(ids));
                  id_winner(j) = ids(id_winner(j));
                else
                  class_scores{j} = [];
                end
              end
            end
            
            [ap(i),fp,tp] = voc_cls_ap(obj.categories{i}, class_scores, obj.imgset);
            dont_care = cellfun(@isempty, class_scores);
            
            n_fp = sum(fp(~dont_care));
            n_tp = sum(tp(~dont_care));
            class_fp{i} = zeros(n_fp,1);
            class_tp{i} = zeros(n_tp,1);
            
            fp_counter = 1;
            tp_counter = 1;
            for j=1:numel(scores)
              if(~dont_care(j))
                if(fp(j))
                  class_fp{i}(fp_counter) = global_ids{j}(id_winner(j));
                  fp_counter = fp_counter + 1;
                elseif (tp(j))
                  class_tp{i}(tp_counter) = global_ids{j}(id_winner(j));
                  tp_counter = tp_counter + 1;
                end
              end
            end                         
          end          
        end      
        
        function [covering_ranking_score, overlap_ranking_score, categ_ranking_score] = eval_ranking_quality(obj, local_ids, MAX_N)
            DefaultVal('*MAX_N', 'inf');
            
            if(MAX_N~=inf)
                for i=1:numel(local_ids)
                    if(~isempty(local_ids{i}))
                        local_ids{i} = local_ids{i}(1:min(numel(local_ids{i}), MAX_N));
                    end
                end
            end
            
            Q = obj.img_overlaps;
            for i=1:numel(local_ids)
                for j=1:numel(Q{i});
                    Q{i}(j).q = Q{i}(j).q(local_ids{i});
                    Q{i}(j).sz = Q{i}(j).sz(local_ids{i});
                end
            end
            
            % the next function is prepared for handling multiple GT
            % segmentations. Must convert to cell.
            for i=1:numel(Q)
                Q{i} = {Q{i}};
            end
                        
            [covering_ranking_score, overlap_ranking_score] = study_segment_quality_nofiles(Q);
            
            % assume there are 20 categories
            for i=1:20
                categ_ranking_score(i) = study_segment_quality_nofiles(Q, i);
            end            
        end         
        
        function Imgs = voc_segm_outputs(obj, img_ids, segm_ids, labels, scores, name, nofiles, algo)                    
            % paste detected segments on the image
            
            DefaultVal('*nofiles', 'false');
            DefaultVal('*algo', '''simple''');
            
            %dir = [obj.exp_dir '/results/' name '/'];
            %if(~exist(dir))
            %    mkdir(dir);
            %end
            
            if(nofiles)
                Imgs = cell(numel(img_ids),1);
            end
            
            cmap = VOClabelcolormap();            
                                        
            parfor i=1:numel(img_ids)
                % load masks
                masks = myload([obj.exp_dir '/MySegmentsMat/' obj.mask_type '/' obj.img_names{img_ids(i)} '.mat'], 'masks');  

                local_segm_ids = obj.global_2_local_whole_ids(segm_ids{i});
                
                newI = zeros(size(masks,1), size(masks,2));                
                [val, ids] = sort(scores{i}, 'ascend');
                
                if(size(ids,1)<=1)
                    sel_masks = masks(:,:,local_segm_ids(ids));
                else
                    % get all (segment combination mode)
                    sel_masks = masks(:,:,local_segm_ids);
                end

                if(strcmp(algo, 'simple'))                
                    for j=1:numel(ids)
                        newI(masks(:,:,local_segm_ids(ids(j)))) = labels{i}(ids(j));
                    end                                    
                else
                    error('no such segment2pixel conversion algorithm');
                end

                if(nofiles)                    
                    Imgs{i} = newI+1;
                else                                                       
                    resfile = sprintf(obj.VOCopts.seg.clsrespath,name,obj.VOCopts.testset,obj.img_names{img_ids(i)});
                    dir_name = resfile(1:end-numel(obj.img_names{img_ids(i)})-numel('.png'));
                    if(~exist(dir_name, 'dir'))
                        mkdir(dir_name);
                    end
                    imwrite(newI+1, cmap, resfile);
                end
            end
        end
        
        function voc_det_outputs(obj, img_ids, segm_ids, labels, scores, name) 
            bndboxes = obj.get_bboxes_img_ids(img_ids);

            the_dir = [obj.exp_dir 'results/VOC2012/Main/'];
            if(~exist(the_dir, 'dir'))
                mkdir(the_dir);
            end            
            
            for i=1:numel(obj.categories)-1 % last one is background
                filename = sprintf([the_dir '%s_det_val_%s.txt'], name, obj.categories{i});
                f = fopen(filename, 'w');                
                for j=1:numel(img_ids)   
                    this_label = find(i == labels{j});
                    ids_this_class = obj.global_2_local_whole_ids(segm_ids{j}(i == labels{j}));
                    for k=1:numel(ids_this_class)
                        m = ids_this_class(k);
                        fprintf(f, '%s %f %f %f %f %f\n', obj.img_names{img_ids(j)}, scores{j}(this_label(k)), bndboxes{j}(m,1), bndboxes{j}(m,2), bndboxes{j}(m,3), bndboxes{j}(m,4));
                    end
                end
                fclose(f);
            end
            disp('ola');
        end

        % Segment quality
        
        function study_segment_quality(obj)
          SvmSegm_study_segment_quality(obj.exp_dir, obj.mask_type, obj.imgset, obj.overlap_type)
        end
    end
end
