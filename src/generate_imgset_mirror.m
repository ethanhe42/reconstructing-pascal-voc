function generate_imgset_mirror(exp_dir, imgset)
    imgset_files = textread([exp_dir 'ImageSets/Segmentation/' imgset '.txt'], '%s');

    f = fopen([exp_dir 'ImageSets/Segmentation/' imgset '_mirror.txt'], 'w');
    for i=1:numel(imgset_files)
        fprintf(f, '%s_mirror\n', imgset_files{i});
    end
    fclose(f);
end