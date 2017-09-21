function SvmSegm_segment_quality_all_JL(img_names, exp_dir, mask_type, segm_quality_type, object_dir, class_dir, break_disc_gt) 
    DefaultVal('*object_dir', '''SegmentationObject''');
    DefaultVal('*class_dir', '''SegmentationClass''');   
    DefaultVal('break_disc_gt', 'false');
    
    segm_eval_dir = [exp_dir 'SegmentEval/'];

    if(~exist(segm_eval_dir))
        mkdir(segm_eval_dir);
    end

    %best_seg = [segm_eval_dir segm_quality_type '/BestSegment/'];
    %if ~exist(best_seg, 'dir')
    %  mkdir(best_seg);
    %end

    if(iscell(mask_type))
        mask_type = cell2mat(mask_type);
    end

    if(break_disc_gt)
        suffix = 'break_disc_gt';
    else
        suffix = '';
    end
    if(~strcmp(object_dir, 'SegmentationObject') || ~strcmp(class_dir, 'SegmentationClass'))
        % non-standard behavior!! create folder !!
        direct = [exp_dir 'SegmentEval/' mask_type '/' segm_quality_type '/' object_dir '_' hifenize(class_dir, suffix) '/'];
        mkdir(direct);
    else        
        direct = [segm_eval_dir mask_type '/' segm_quality_type '/' suffix '/'];
        if ~exist(direct, 'dir')
            mkdir(direct);
        end
    end    
    
    t1 =  tic();
    for i=1:length(img_names)    
        if(exist([direct img_names{i} '.mat'], 'file'))
            continue;
        end

        masks = myload([exp_dir 'MySegmentsMat/' mask_type '/' img_names{i} '.mat'], 'masks');
        if iscell(masks)
          masks = cell2mat(masks);
        end

        name = [exp_dir object_dir '/' img_names{i} '.png'];
        alternatives = dir([name '*']);

        name =  [exp_dir class_dir '/' img_names{i} '.png'];
        class_alternatives = dir([name '*' ]);

        assert(~isempty(alternatives));
        assert(numel(alternatives)==numel(class_alternatives));

        ground_truth_obj_segs = cell(numel(alternatives),1);

        Quality = cell(1, numel(alternatives));
        for j=1:numel(alternatives)
            %Quality = struct();

            % read class image
            [classI, map] = imread([exp_dir class_dir '/' class_alternatives(j).name]);

            % read instance image 
            ground_truth_obj_segs{j} = imread([exp_dir object_dir '/'  alternatives(j).name]);

            if(break_disc_gt)
                labels= unique(ground_truth_obj_segs{j});

                next_label = max(labels)+1;
                for k =1:numel(labels)
                    p = bwconncomp(ground_truth_obj_segs{j}== k);
                    if(numel(p.PixelIdxList)>1)
                        for l=2:(numel(p.PixelIdxList))
                            ground_truth_obj_segs{j}(p.PixelIdxList{l}) = next_label;
                            next_label = next_label +1;
                        end
                    end
                end
            end
            
            un = unique(ground_truth_obj_segs{j});
            un(un==0) = [];
            un(un==255) = [];
            care = (classI~=255);

            for k=1:numel(un)                
                ground_truth_k = zeros(size(ground_truth_obj_segs{j}));
                ground_truth_k(ground_truth_obj_segs{j} == un(k)) = 1;

                %disp('debugginnnnnn')
                %sz_masks = size(masks)
                %sz_gt = size(ground_truth_k)
                [duh1, duh2, qual] = myCalcCandScoreFigureGroundAll(masks,ground_truth_k, segm_quality_type, care);

                new_un = unique(ground_truth_k);
                new_un(new_un==0) = [];
                count = [];
                [clslbl, a,c_occ] = unique(classI(ground_truth_k==new_un)); % can do faster than this
                if(numel(clslbl)~=1) % noisy dataset, single object region has multiple classes
                    disp('noisy dataset!');
                    % get the lbl that appear most
                    count = [];
                    for l=1:numel(clslbl)
                        count(l) = sum(c_occ==l);
                    end
                    [a, id_max] =max(count);
                    clslbl = clslbl(id_max);
                end

                %img_names{i}
                %clslbl
                %return;
                Quality{j}(k).class = clslbl;
                Quality{j}(k).object = k;
                Quality{j}(k).object_sz = sum(sum(ground_truth_k));
                Quality{j}(k).image_name = img_names{i};

                if(isfield(Quality{j}(k), 'q') && ~isempty(Quality{j}(k).q))
                    Quality{j}(k).q
                    error('test this first');
                    for l=1:numel(Quality{j}(k).q)
                        if(Quality{j}(k).q(l) < qual(l))
                            Quality{j}(k).q(l) = qual(l);
                            Quality{j}(k).sz(l) = sum(sum(masks(:,:,l)));
                        end
                    end
                    Quality{j}(k).sz = sum(masks,3);
                else
                    Quality{j}(k).q = qual;
                    Quality{j}(k).sz = squeeze(sum(sum(masks,1), 2));
                end
            end
        end

        file_to_save = [direct img_names{i} '.mat'];
        mysave(file_to_save, 'Quality', Quality);
        %save(file_to_save, 'Quality');
    end
end
