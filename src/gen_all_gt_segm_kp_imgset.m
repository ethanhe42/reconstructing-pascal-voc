% generate imgset containing only those images where all objects have
% keypoints and segmentations and where there's no obvious problem 
% with the keypoints (eg. keypoints missing on one object, repeated for another)
function filename = gen_all_gt_segm_kp_imgset()    
    exp_dir = add_all_paths();
    mask_type = 'ground_truth';

    base_imgset = 'all_gt_segm';

    sb = SegmBrowser(exp_dir, mask_type, base_imgset);
    corr = CorrespBRKL(exp_dir, base_imgset, sb);

    n_with_annot = zeros(numel(sb.img_names),1);
    n_obj = zeros(numel(sb.img_names),1);
    n_gt_segm = n_obj;
    selected = false(numel(sb.img_names),1);
    for i=1:numel(sb.img_names)
        i

        %I = sb.get_Imgs(i);
        %I = I{1};

        count_with_kp = arrayfun(@(a) numel(a.keypoints), corr.img_ann(i).objects);
        n_with_annot(i) = sum(count_with_kp ~= 0);
        n_obj(i) = numel(count_with_kp);

        n_gt_segm(i) = numel(sb.img_2_whole_ids{i});
        %[F, D] = pooling_local_feats_mask_sift_kp(I, masks, kp, pars)

        selected(i) = false;
        if sum(n_with_annot(i)==n_obj(i))
           % do more expensive test

           if(1)
               masks = sb.get_img_masks(i);
                               
               this_dist = zeros(numel(corr.img_ann(i).objects),1);
               for j=1:numel(corr.img_ann(i).objects)
                   kp = ceil(corr.img_ann(i).objects(j).keypoints);
                   visib = logical(corr.img_ann(i).objects(j).visib);
                   kp_ids = find(kp(:,1)~=0.0 | kp(:,2)~=0.0);
                   
                   assert(numel(corr.img_ann(i).objects(j).visib) == numel(kp_ids));
                   
                   kp_visib = kp(kp_ids(visib),:);
                   kp_visib(kp_visib(:,1) == 0 | kp_visib(:,2) == 0, :) = [];
                   
                   dist = bwdist(masks(:,:,j));
                   
                   try
                      %MERCIFULNESS = 0;
                      %kp_visib(kp_visib<0 & kp_visib >-MERCIFULNESS) = 1;
                      %kp_visib(kp_visib(:,1) > (size(masks,2)+MERCIFULNESS),1) = size(masks,2);
                      %kp_visib(kp_visib(:,2) > (size(masks,1)+MERCIFULNESS),1) = size(masks,1);
                      dists = dist(sub2ind([size(masks,1) size(masks,2)], kp_visib(:,2), kp_visib(:,1)));
                      this_dist(j) = mean(dists);
                   catch                       
                       this_dist(j) = inf;
                   end
               end

               THRESH = 15;
               if(max(this_dist)<THRESH)
                   selected(i) = true;
               end

               if(~selected(i) && 0)
                   % visualize
                   close all;

                   corr.show_wholes_and_keypoints(sb.img_2_whole_ids{i});
                   disp('oh oh');
                   figure;

                   I = imread([exp_dir 'JPEGImages/' sb.img_names{i} '.jpg']);
                   imshow(I); hold on;
                   for j=1:numel(corr.img_ann(i).objects)
                      kp = corr.img_ann(i).objects(j).keypoints;
                      kp(kp(:,1) == 0 & kp(:,2) == 0,:) = [];
                      plot(kp(:,1), kp(:,2), 'o', 'LineWidth', 15, 'Color', cmap(j+1,:));
                      axis image;
                   end
               end
           end
        end
    end  

    missing_ann_objs = sum(n_obj) - sum(n_with_annot)

    %ratio_without_kp = sum(n_no_annot)/sum(n_obj)
    ratio_zero_kp_img = sum(n_with_annot==n_obj)/numel(n_obj)
    sel_imgs = sb.img_names(selected);
    assert(all(n_obj(n_with_annot==n_obj) == n_gt_segm(n_with_annot==n_obj)))

    if 1
        filename = [base_imgset '_kp.txt'];
        f = fopen([exp_dir 'ImageSets/Segmentation/' filename], 'w');
        for i=1:numel(sel_imgs)
            fprintf(f, '%s\n', sel_imgs{i});
        end
        fclose(f);
    end
