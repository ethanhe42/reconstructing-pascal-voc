% Joao Carreira September 2012
% Run this in debug mode and do step by step for checking the outputs it gives

function test_SegmBrowser()
  % run this after setting up the VOC dataset (see in VOC_experiment
  % folder)
  exp_dir = '../../VOC_experiment/VOC/';
  mask_type = 'CPMC_segms_150_sp_approx';
  imgset = 'train11';  
  overlap_type = 'overlap';  
  
  addpath('../');
  addpath('../../external_src/');
  addpath('../../external_src/immerge/');
  
  % first time it runs it may take a few seconds, then it caches everything
  % it needs
  b = SegmBrowser(exp_dir, mask_type, imgset, overlap_type);
  
  %
  % Images, objects and segments have numerical identifiers that apply
  % given the imgset
  % 
  % If you want to get features for both ThreeAmigos and GroundTruth masks,
  % I suggest you create two SegmBrowser objects, one for each type.
  %
    
  imgname = '2007_000032';
  
  % get id of the image 
  img_id = b.img_names_to_ids(imgname);
    
  % show image
  b.show_imgs(img_id);  
  
  if 1
    % get segment (whole) ids for that image
    segment_ids = b.img_2_whole_ids(img_id);    
    % show first 50 segments for image
    b.show_wholes(segment_ids{1}(1:50))    
    
    % show best segments for that image
    b.show_best_masks(imgname);
    
    % get overlaps for image
    overlap_img = b.get_overlaps_img(b.img_names{img_id});      
  end
  
  if 1
    % get labels for all segments (by finding max overlap)
    [labels] = b.get_labels_wholes();
    % get img object ids and overlaps for a particular class, one cell for
    % each image
    [img_ids, obj_ids, Q] = b.collect_category_imgs(1) 
    % show objects in first image of this category
    b.show_objects(obj_ids{1});  
    % get overlap for all segments of class 1
    [overlap_class] = b.get_overlap_wholes_objclass(1);
  end
  
  % to get all best masks for cow, not sure if useful
  categ_masks = b.collect_all_category_best_masks('cow');     
end