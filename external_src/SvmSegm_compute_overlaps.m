function SvmSegm_compute_overlaps(exp_dir, img_names, mask_type, force)
    if(~exist('force', 'var'))
        force = false;
    end

    overlaps_dir = [exp_dir 'MyOverlaps/'];
    masks_dir = [exp_dir 'MySegmentsMat/'];

    if(~exist([overlaps_dir mask_type], 'dir'))
        mkdir([overlaps_dir mask_type]);
    end

    for i=1:length(img_names)
        out_name = [overlaps_dir mask_type '/' img_names{i} '.mat'];
        if(exist(out_name, 'file') && ~force)
            continue
        end
               
        masks = logical([]);        
 
        filename = [masks_dir mask_type '/' img_names{i} '.mat'];
        if(~exist(filename, 'file'))
            filename
            error('masks missing');
        else
            var = load(filename);
            masks = cat(3, masks, var.masks);
            if(isempty(masks))
                error('no masks inside!');
            end
        end
        
        
        masks = imresize(masks, 0.25, 'nearest');
        masks = reshape(masks, size(masks,1) * size(masks,2), size(masks,3));
        overlap_mat = single(segm_overlap_mex(masks));
        save(out_name, ['overlap_mat']);
    end
end
