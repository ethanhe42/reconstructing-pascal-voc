function download_pascal_dataset(exp_dir)    
    % obtain the most recent pascal voc annotations
    voc_file = 'VOCtrainval_11-May-2012.tar';
    voc_link = 'http://pascallin.ecs.soton.ac.uk/challenges/VOC/voc2012/VOCtrainval_11-May-2012.tar';

    mkdir(exp_dir);
    mkdir([exp_dir 'SegmentEval/']);

    % Desired path to VOC images
    img_path = [exp_dir 'JPEGImages/'];
    % We will train using ground truth annotations for all images, with
    % the additional Berkeley annotations

    if(~exist('./benchmark.tgz', 'file'))
        % Get the Berkeley annotations
        disp('Downloading external voc ground truth annotations.');
        disp('See http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/semantic_contours/');
        pause(5);
        
        !wget http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/semantic_contours/benchmark.tgz

        % Unpack Berkeley annotations, store in VOC format
        !tar xfz benchmark.tgz

        files = dir('./benchmark_RELEASE/dataset/inst/*.mat');
        mkdir([exp_dir 'SegmentationObject/']);
        mkdir([exp_dir 'SegmentationClass/']);
        cmap = VOClabelcolormap(256);
        parfor i=1:numel(files)
            var = load(['./benchmark_RELEASE/dataset/inst/' files(i).name]);
            I_obj = var.GTinst.Segmentation;
            I_cls = I_obj;
            for j=1:numel(var.GTinst.Categories)
                I_cls(I_obj==j) = var.GTinst.Categories(j);
            end

            I_obj = I_obj+1;
            I_cls = I_cls+1;
            imwrite(I_obj, cmap, [exp_dir 'SegmentationObject/' files(i).name(1:end-4) '.png']);
            imwrite(I_cls, cmap, [exp_dir 'SegmentationClass/' files(i).name(1:end-4) '.png']);
        end
    end

    if(~exist(voc_file, 'file'))
        % Get the latest higher-quality Pascal VOC annotations and replace
        % any overlapping Berkeley annotations by these
        disp(['Downloading PASCAL VOC dataset. See  http://pascallin.ecs.soton.ac.uk/challenges/VOC/']);
        pause(5);
        
        system(['wget ' voc_link]);
        system(['tar xf ' voc_file]);
        mkdir(img_path);
        voc_dir = dir('VOCdevkit');
        voc_dir(1:2) = [];
        system(['cp VOCdevkit/' voc_dir.name '/JPEGImages/* ' img_path]);
        system(['cp VOCdevkit/' voc_dir.name '/SegmentationObject/* ' exp_dir 'SegmentationObject/']);
        system(['cp VOCdevkit/' voc_dir.name '/SegmentationClass/* ' exp_dir 'SegmentationClass/']);
        system(['cp -r VOCdevkit/' voc_dir.name '/Annotations/ ' exp_dir]);
    end

    % Copy imgset files
    voc_dir = dir('VOCdevkit');
    voc_dir(1:2) = [];
    system(['cp -r VOCdevkit/' voc_dir.name '/ImageSets/ ' exp_dir]);

    % Create all_gt_segm imgset and variations
    imgset_dir = [exp_dir 'ImageSets/Segmentation/'];
    list_all_gt = dir([exp_dir 'SegmentationObject/*.png']);
    list_all_gt = {list_all_gt(:).name};
    list_all_gt_segm = cellfun(@(a) a(1:end-4), list_all_gt, 'UniformOutput', false);
    f = fopen([imgset_dir 'all_gt_segm.txt'], 'w');
    for i=1:numel(list_all_gt_segm)
        fprintf(f, '%s\n', list_all_gt_segm{i});
    end
    fclose(f);
    generate_setdiff_imgset(exp_dir, 'all_gt_segm', 'val');
    generate_imgset_mirror(exp_dir, 'all_gt_segm');
    generate_imgset_mirror(exp_dir, 'all_gt_segm_minus_val');
    generate_imgset_mirror(exp_dir, 'train');
    generate_imgset_mirror(exp_dir, 'val');
    generate_imgset_mirror(exp_dir, 'trainval');

    % Generate ground truth masks for all images
    disp('Generating ground truth masks.');    
    img_names = textread([exp_dir 'ImageSets/Segmentation/all_gt_segm.txt'], '%s');
    SvmSegm_generate_gt_masks(exp_dir, img_names);

    disp('Creating mirrored images.');
    create_mirrored_images(exp_dir, 'all_gt_segm')
    mirror_imgset = 'all_gt_segm_mirror';

    disp('Generating ground truth masks for mirrored images.');
    img_names = textread([exp_dir 'ImageSets/Segmentation/' mirror_imgset '.txt'], '%s');
    SvmSegm_generate_gt_masks(exp_dir, img_names);

    % Compute all necessary segment metadata
    a = SegmBrowser(exp_dir, 'ground_truth', 'all_gt_segm');
    b = SegmBrowser(exp_dir, 'ground_truth', 'all_gt_segm_mirror');
    
    
    %%%%%%%%%%%% Get keypoint annotations %%%%%%%%%%%
    disp('Downloading external voc keypoint annotations.');
    disp('See http://www.cs.berkeley.edu/~lbourdev/poselets/');
    !wget http://www.cs.berkeley.edu/~lbourdev/poselets/voc2011_keypoints_Feb2012.tgz
    % goes to folder "annotations"
    !tar xfz voc2011_keypoints_Feb2012.tgz 
    
    get_kp_names();
    convert_dataset_kp_format();
    collect_imgset_keypoints(exp_dir, {'all_gt_segm'});
    
    system(['cp ./annotations/unique* ' exp_dir '/Browser/']);
    assemble_voc_kp_metadata(); % gather left-right symmetry metadata
    collect_imgset_keypoints_mirror(exp_dir, {'all_gt_segm'}, {'all_gt_segm_mirror'});
    
    % generate imgset with those images having clean annotations for all
    % pascal voc objects in them
    gen_all_gt_segm_kp_imgset();
    generate_imgset_mirror(exp_dir, 'all_gt_segm_kp');
    system(['mv ' exp_dir 'ImageSets/Segmentation/all_gt_segm_kp_mirror.txt ' exp_dir 'ImageSets/Segmentation/all_gt_segm_mirror_kp.txt']);
    
    % Compute segment metadata for these image sets
    a = SegmBrowser(exp_dir, 'ground_truth', 'all_gt_segm_kp');
    b = SegmBrowser(exp_dir, 'ground_truth', 'all_gt_segm_mirror_kp');
    
    collect_imgset_keypoints(exp_dir, {'all_gt_segm_kp'});
    collect_imgset_keypoints_mirror(exp_dir, {'all_gt_segm_kp'}, {'all_gt_segm_mirror_kp'});
end

