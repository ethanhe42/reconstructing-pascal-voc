function SvmSegm_generate_gt_masks(exp_dir, img_names, filter_small)
    DefaultVal('*filter_small', 'false');
    gt_folder = [exp_dir 'SegmentationObject/'];
    gt_class_folder = [exp_dir 'SegmentationClass/'];
    segm_folder = [exp_dir 'MySegmentsMat/'];    
    
    SMALL = 50; % less than this it gets filtered out, if filter small is on
    
    counter = 0;
    
    %img_names = textread([exp_dir 'ImageSets/Segmentation/' imgset '.txt'], '%s');
     
    for i=1:numel(img_names)
        if(exist([exp_dir 'MySegmentsMat/ground_truth/' img_names{i}  '.mat']))
            continue;
        end
                
        % disable this functionality ( it is useful for bsds )
        %alternatives = dir([gt_class_folder img_names{i} '*']); 
        alternatives = dir([gt_class_folder img_names{i} '.png']);        
        I = imread([gt_folder alternatives(1).name]);
        if(~exist([segm_folder '/ground_truth/'], 'dir'))
            mkdir([segm_folder '/ground_truth/']);
        end
        %masks = cell(numel(alternatives), 1);
        counter = 1;
        
        masks = false(size(I,1), size(I,2), 200);
        gt_tiling_ids = [];
        for j=1:numel(alternatives)
            I = imread([gt_folder alternatives(j).name]);
            un = unique(I);
            un(un==0) = [];
            un(un==255) = []; % in voc 255 is for hard or don't care
            
            %masks{j} = false(size(I));

            for k=1:numel(un)
                newI = false(size(I));
                newI(I==un(k)) = true;
                if(filter_small && (SMALL>sum(newI(:))))
                    continue;
                end
                masks(:,:,counter) = newI;
                gt_tiling_ids(counter) = j;
                counter = counter + 1;
                %if(k==1)
                    %masks{j} = newI;                  
                %else
                %    masks(:,:,j) = 
                    %masks{j} = cat(3, masks{j}, newI);
                %end
            end
        end            
        masks(:,:,counter:end) = [];
        save([exp_dir 'MySegmentsMat/ground_truth/' img_names{i}  '.mat'], 'masks', 'gt_tiling_ids');
    end    
end