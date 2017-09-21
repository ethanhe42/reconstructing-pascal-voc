function SvmSegm_segment_quality_JL_join(exp_dir, data_partitions, mask_type, segm_quality_type, another_challenge, segm_gt_type)
  DefaultVal('*segm_gt_type', '[]');  
  DefaultVal('*another_challenge', '[]');
  segm_eval_dir = [exp_dir 'SegmentEval/'];

  direct = [segm_eval_dir mask_type '/' segm_quality_type '/' segm_gt_type '/'];

  for h=1:numel(data_partitions)
    %img_names =textread(sprintf(opts.seg.imgsetpath,data_partitions{h}), '%s');
    if exist('another_challenge') && ~isempty(another_challenge);
        img_names = textread([exp_dir 'ImageSets/' another_challenge '/' data_partitions{h} '.txt'], '%s');
    else
        img_names = textread([exp_dir 'ImageSets/Segmentation/' data_partitions{h} '.txt'], '%s');
    end
    %Q = zeros(numel(img_names),1);
    Q = [];
    for i=1:numel(img_names)
        try
          load([direct img_names{i}]);        
        catch
          if(~strcmp(segm_quality_type, 'overlap_withBG'))
            SvmSegm_segment_quality_all_JL(img_names, exp_dir, mask_type, segm_quality_type, another_challenge, segm_gt_type);
          else
            SvmSegm_overlap_with_background(img_names(i), exp_dir, mask_type, 'overlap');
          end
          
          load([direct img_names{i}]);
        end
        
        if(~iscell(Quality));
            Quality = {Quality};
        end
        
        newQuality = []; % merge different ground truth segmentations
        for j=1:numel(Quality)
            newQuality = [newQuality Quality{j}];
        end
        Q = [Q newQuality]; % always stores the first (might need to have a look here)
    end
    Quality = Q;
    save([direct data_partitions{h} '.mat'], 'Quality', '-V6'); % -V6 -- no compression (10x faster!)
  end
end
