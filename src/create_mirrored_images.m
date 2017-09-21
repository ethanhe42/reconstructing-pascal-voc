function create_mirrored_images(exp_dir, imgset)    
    imgset_dir = [exp_dir 'ImageSets/Segmentation/'];
    img_names = textread([imgset_dir imgset '.txt'], '%s');
    
    % generate any missing .jpg files
    JPEGdir = [exp_dir 'JPEGImages/'];
    for i=1:numel(img_names)
        jpeg_name = [JPEGdir img_names{i} '_mirror.jpg'];
        if(~exist([jpeg_name], 'file'))
            % mirror the image
            I = imread([JPEGdir img_names{i} '.jpg']);
            I = I(:, abs(-size(I,2):-1), :);            
            imwrite(I, jpeg_name);
        end
    end
    
    % generate respective object and class annotations
    cls_dir = [exp_dir 'SegmentationClass/'];
    obj_dir = [exp_dir 'SegmentationObject/'];
    for i=1:numel(img_names)
        png_cls_name = [cls_dir img_names{i} '_mirror.png'];
        png_obj_name = [obj_dir img_names{i} '_mirror.png'];        
        
        % mirror the image
        [Icls,cmap] = imread([cls_dir img_names{i} '.png']);
        Icls = Icls(:, abs(-size(Icls,2):-1), :);
        imwrite(Icls, cmap, png_cls_name);
        
        [Iobj, cmap] = imread([obj_dir img_names{i} '.png']);
        Iobj = Iobj(:, abs(-size(Iobj,2):-1), :);
        imwrite(Iobj, cmap, png_obj_name);        
    end
end