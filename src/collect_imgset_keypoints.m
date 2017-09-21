function collect_imgset_keypoints(exp_dir, imgset)
    dest_dir = [exp_dir 'merged_Correspondences_GT_BRKL/'];
    
    if(~exist(dest_dir, 'dir'))
        mkdir(dest_dir);
    end
    
    total_n_objects = 0 ;
        
    for h=1:numel(imgset)
        img_names = textread([exp_dir 'ImageSets/Segmentation/' imgset{h} '.txt'], '%s');    
        
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
            
            try 
                GT_rec_ours = load([exp_dir 'SegmentEval/ground_truth/overlap/' img_names{i} '.mat']);
                GT_rec_ours = GT_rec_ours.Quality{1};
                if(numel(var.rec.objects) ~= numel(GT_rec_ours))
                    % here we fix the dessynchronization between PASrecords and our Quality files (see your research notes)!!!
                    var.rec.objects = var.rec.objects(1:numel(GT_rec_ours));                    
                end
            catch
                disp('no such gt file');
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
        
        save([dest_dir imgset{h} '.mat'], 'rec', 'obj_keypoints', '-V6');
    end
end
 