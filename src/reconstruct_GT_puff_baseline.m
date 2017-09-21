function reconstruct_GT_puff_baseline(exp_dir, name, sel_class)        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Get data of images having ground truth meshes %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    imgset_all = 'synth_5_views_per_mesh';
    synth_mask = 'ground_truth';
    sb_synth = SegmBrowser(exp_dir, synth_mask, imgset_all);
    y = sb_synth.get_overlaps_wholes(1:numel(sb_synth.whole_2_img_ids));
    [q, class] = max(y, [], 2);
    in_class = find(class == sel_class);    
    masks = sb_synth.get_masks(in_class);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Reconstruct object instances %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %all_reconstr = repmat(struct('vertices', single(0), 'faces', int32(0),'dt', single(0)), size(in_class,1), 1);
    t_total = tic();
    
    out_folder = ['Results/' name '/' sb_synth.categories{sel_class} '/'];
    mkdir(out_folder);
    
    save([out_folder 'reconstruction_data'], 'in_class', 'imgset_all','-v7.3');
    parfor i=1:numel(masks)        
        tic
        disp(i);
        file_name = [out_folder int2str(i) '.mat'];
        if(exist(file_name,'file'));
           continue
        end

        reconstructions = [];
        [reconstructions.faces, reconstructions.vertices] = puffball(masks{i});   
        reconstructions.vertices = reconstructions.vertices';
        reconstructions.score = 1;
        reconstructions.triples = i;

        if 0
            I = sb_synth.get_Imgs(sb_synth.whole_2_img_ids(in_class));
            I = I{1};
            
            % visualize reconstruction 
            sc(I); hold on;
            
            trisurf(reconstructions.faces, reconstructions.vertices(1,:), ...
                reconstructions.vertices(2,:), reconstructions.vertices(3,:), 'FaceColor', 'red');
            axis equal;
            set(gca,'zdir','reverse');
            close all;
        end
        mysave(file_name, 'reconstructions', reconstructions);
        toc
    end
    time = toc(t_total)
    save([out_folder 'reconstruction_data'], 'in_class', 'imgset_all', 'time','-v7.3');
end
