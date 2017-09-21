function [] = reconstruct_GT_convexHull_baseline(sel_class, exp_dir,n_iter_1,n_iter_2,name)

imgset_GT = 'synth_5_views_per_mesh';
mask_type = 'ground_truth';
sb_synth = SegmBrowser(exp_dir, mask_type, imgset_GT);

imgset_all = 'all_gt_segm_kp_synth_5_views_per_mesh';
sb_all = SegmBrowser(exp_dir, mask_type, imgset_all);


load voc_kp_metadata

sfm_folder = ['./SFM_models/'];

if(metadata.articulated(sel_class))
    largest_rigid = true;        
else
    largest_rigid = false;
end


[in_class, filename_sfm, filename_ref] = compute_and_cache_sfm_model(exp_dir, sfm_folder, sel_class, n_iter_1, n_iter_2, largest_rigid, imgset_all);
load([filename_sfm '.mat'],'Shape');
load([filename_ref '.mat'], 'R', 'T');        
[Shape, R] = postprocess_sfm(Shape, R);



y = sb_synth.get_overlaps_wholes(1:numel(sb_synth.whole_2_img_ids));
[q, class] = max(y, [], 2);
in_class_synth = find(class == sel_class);


ids_all = sb_all.whole_2_img_ids(in_class);
[~,synth_ids,~] = intersect(sb_all.img_names(ids_all), sb_synth.img_names(in_class_synth));



tri_CH = convhull(Shape(1,:)',Shape(2,:)',Shape(3,:)');

out_folder = ['Results/' name '/' metadata.categories{sel_class} '/'];
mkdir(out_folder);


save([out_folder 'reconstruction_data'], 'in_class', 'imgset_all');

for i=1:numel(synth_ids)        
    disp(i);
    file_name = [out_folder int2str(i) '.mat'];

    theR = R(:,:,synth_ids(i));
    theT = T(:,synth_ids(i));
    
    reconstructions = [];
    reconstructions.faces = tri_CH;
    reconstructions.vertices = bsxfun(@plus,theR*Shape,theT);
    reconstructions.score = 1;
    reconstructions.triples = synth_ids(i);    

    save(file_name, 'reconstructions');
end


end

