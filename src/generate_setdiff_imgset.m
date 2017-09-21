function generate_setdiff_imgset(exp_dir, imgset1, imgset2)    
    imgset_out = [exp_dir 'ImageSets/Segmentation/' imgset1 '_minus_' imgset2];
        
    files_1 = textread([exp_dir 'ImageSets/Segmentation/' imgset1 '.txt'], '%s');
    files_2 = textread([exp_dir 'ImageSets/Segmentation/' imgset2 '.txt'], '%s');
    
    files_out = setdiff(files_1, files_2);
        
    filename = [imgset_out '.txt'];
    f = fopen(['./' filename], 'w');
    for i=1:numel(files_out )
        fprintf(f, '%s\n', files_out{i});
    end
    fclose(f);    
end