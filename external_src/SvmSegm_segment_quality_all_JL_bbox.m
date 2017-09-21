function SvmSegm_segment_quality_all_JL_bbox(img_names, exp_dir, mask_type, segm_quality_type)
assert(strcmp(segm_quality_type, 'overlap_and_bbox_edges'));
VOCinit();

segm_eval_dir = [exp_dir 'SegmentEval/'];

if(~exist(segm_eval_dir))
    mkdir(segm_eval_dir);
end

%best_seg = [segm_eval_dir segm_quality_type '/BestSegment/'];
%if ~exist(best_seg, 'dir')
%  mkdir(best_seg);
%end

direct = [segm_eval_dir mask_type '/' segm_quality_type '/' ];
if ~exist(direct, 'dir')
    mkdir(direct);
end

t1 =  tic();
for i=1:length(img_names)
    %for i=1:length(img_names)
    if(exist([direct img_names{i} '.mat'], 'file'))
        continue;
    end
    
    masks = myload([exp_dir 'MySegmentsMat/' mask_type '/' img_names{i} '.mat'], 'masks');
    if iscell(masks)
        masks = cell2mat(masks);
    end
    
    ground_truth_obj_segs = cell(1,1);
    
    Quality = cell(1, 1);
    alternatives = 1;
    for j=1:numel(alternatives)
        rec=PASreadrecord(sprintf(VOCopts.annopath,img_names{i}));
       
        %to_remove = false(numel(rec.objects),1);
        %for k=1:numel(rec.objects)
        %    if(rec.objects(k).difficult)
        %        to_remove(k) = true;
        %    end
        %end
        %rec.objects(to_remove) = [];
        
        n_objects = numel(rec.objects);
        if(n_objects ==0)
            % add dummy object, set all qualities to zero
            n_objects = 1;
            rec.objects(1).bbox = [-10 -10 -5 -5];
            rec.objects(1).class = {'dummy'};
        end
        
        care = ones(size(masks,1), size(masks,2));
        for k=1:n_objects
            k
            ground_truth_k = false(size(masks,1), size(masks,2));
            bbox = rec.objects(k).bbox;
            
            if(rec.objects(1).bbox(1)>0) % not dummy object
                ground_truth_k(bbox(2):bbox(4), bbox(1):bbox(3)) = true;                                
            end
            
            %disp('debugginnnnnn')
            %sz_masks = size(masks)
            %sz_gt = size(ground_truth_k)
            [duh1, duh2, qual] = myCalcCandScoreFigureGroundAll(masks,ground_truth_k, segm_quality_type, care);
            [a, clslbl] = intersect(VOCopts.classes, rec.objects(k).class);
            if(numel(clslbl) ~= 1)
                clslbl = -1;
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
end
toc(t1)

end
