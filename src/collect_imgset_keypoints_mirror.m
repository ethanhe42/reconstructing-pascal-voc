function collect_imgset_keypoints_mirror(exp_dir, imgset, imgset_mirror)
    dest_dir = [exp_dir 'merged_Correspondences_GT_BRKL/'];
    
    load('./voc_kp_metadata.mat', 'metadata');
    
    if(~exist(dest_dir, 'dir'))
        mkdir(dest_dir);
    end
    
    total_n_objects = 0 ;
    
    assert(iscell(imgset));
    assert(numel(imgset) == numel(imgset_mirror));
    
    for h=1:numel(imgset_mirror)
        
        sb = SegmBrowser(exp_dir, 'ground_truth', imgset{h});
        y = sb.get_overlaps_wholes(1:numel(sb.whole_2_img_ids));
        [q, classes] = max(y, [], 2);
        classes = mat2cell(classes, cellfun(@numel, sb.img_2_whole_ids));
        
        img_names = sb.img_names;

        counter = 1;
        obj_keypoints = {};
        all_rec = [];
        for i=1:numel(img_names)                                        
            i
            var = load([dest_dir img_names{i} '.mat'], 'rec');            
            
            if 0
                % debugging
                VOCinit();
                GT_rec = PASreadrecord(sprintf(VOCopts.annopath,img_names{i}));
                GT_rec_ours = load([exp_dir 'SegmentEval/ground_truth/overlap/' img_names{i} '.mat']);
                GT_rec_ours = GT_rec_ours.Quality{1};
                assert(numel(var.rec.objects) == numel(GT_rec_ours));
                total_n_objects = total_n_objects + numel(GT_rec_ours);
            end                        
            
            GT_rec_ours = load([exp_dir 'SegmentEval/ground_truth/overlap/' img_names{i} '.mat']);
            GT_rec_ours = GT_rec_ours.Quality{1};
            if(numel(var.rec.objects) ~= numel(GT_rec_ours))
                % here we fix the dessynchronization between PASrecords and our Quality files (see your research notes)!!!
                var.rec.objects = var.rec.objects(1:numel(GT_rec_ours));                    
            end
            
            if 0 % debug
                I = imread([exp_dir 'JPEGImages/' img_names{i} '.jpg']);
                imshow(I); hold on;
                %Iobj = imread([exp_dir 'SegmentationObject/' img_names{i} '.png']);
                cmap = VOClabelcolormap(256);
                %imshow(Iobj, cmap); hold on;
                for j=1:numel(var.rec.objects)
                    if(isempty(var.rec.objects(j).keypoints))
                        continue;
                    end
                    plot(var.rec.objects(j).keypoints(:,1), var.rec.objects(j).keypoints(:,2), 'o', 'LineWidth', 5, 'Color', cmap(j+1,:));
                    if 0
                        % fix object annotation
                        
                    end
                end                
            end
            
            I = imread([exp_dir 'JPEGImages/' img_names{i} '.jpg']);
            img_size = [size(I,2) size(I,1)];            
            for j=1:numel(var.rec.objects)        
                if(isempty(var.rec.objects(j).keypoints))
                    continue;
                end
                                
                non_null = (sum(var.rec.objects(j).keypoints,2) ~= 0);
                var.rec.objects(j).keypoints(non_null, 1) = img_size(1) - var.rec.objects(j).keypoints(non_null,1) + 1;                
                
                map = metadata.sym_corresp{classes{i}(j)}{1};                                
                var.rec.objects(j).keypoints(:, 1)  = swap(var.rec.objects(j).keypoints(:, 1) , map(1,:), map(2,:));
                var.rec.objects(j).keypoints(:, 2)  = swap(var.rec.objects(j).keypoints(:, 2) , map(1,:), map(2,:));
                % swap left-right symmetric dimensions
            end
            
            all_rec= [all_rec var.rec];
            %all_rec(i) = var.rec;
                        
            for j=1:numel(all_rec(i).objects)
                if(isfield(all_rec(i).objects(j), 'keypoints'))
                    obj_keypoints{counter} = all_rec(i).objects(j).keypoints;
                else
                    obj_keypoints{counter} = [];
                end
                
                counter = counter + 1;
            end
        end
        
        rec = all_rec;
        
        save([dest_dir imgset_mirror{h} '.mat'], 'rec', 'obj_keypoints', '-V6');
    end
end
 
function vector = swap(vector, dims1, dims2)
    vec2 = vector;
    vec2(dims1) = vector(dims2);
    vec2(dims2) = vector(dims1);
    vector = vec2;
end