function get_kp_names()
    exp_dir = add_all_paths();
    mask_type = 'ground_truth';
    imgset = 'all_gt_segm'; % doesn't matter, we'll use all annotations

    warning off;

    sb = SegmBrowser(exp_dir, mask_type, imgset);

    VOCinit();

    kp_names = [];   
    for j=1:numel(sb.categories)
        if(exist(['./annotations/unique_kp_' sb.categories{j} '.mat'], 'file'))
            %continue;
        end
        dir_names = ['./annotations/' sb.categories{j} '/*.xml'];
        filenames = dir(dir_names);
        if(isempty(filenames))
            % can't find files for that class
            disp(['missing files for ' sb.categories{un_classes_img(j)} ' and img ' img_name ' !']);
            continue;
        end
        
        kp_names = cell(numel(filenames),1);
        parfor k=1:numel(filenames)
            %filenames(k).name
            new_struct = read_brkl_kp(['./annotations/' sb.categories{j} '/' filenames(k).name]);
            kp_names{k} = new_struct.kp_name;
        end        
        
        empty = find(cellfun(@isempty, kp_names));
        kp_names(empty) = [];
        all_names = [];
        for k=1:numel(kp_names)
            all_names = [all_names; kp_names{k}];
        end
            
        un_names = unique(all_names);
        save(['./annotations/unique_kp_' sb.categories{j} '.mat'], 'un_names');
    end
    
end