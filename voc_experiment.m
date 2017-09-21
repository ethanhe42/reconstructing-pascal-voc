function voc_experiment()                     
    exp_dir = add_all_paths();      
    eval('config');
    init_rnd_seed();
    
    % compile mex files 
    compile_all();

    % create multiple threads (set how many you have)
    if(matlabpool('size')~=N_THREADS)
        matlabpool('open', N_THREADS);
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    t_download = tic();    
    % 1. Get the data
    if(~exist('VOC/MySegmentsMat/ground_truth/', 'dir')) % do once
        download_pascal_dataset(exp_dir);
    end
    toc(t_download)
    
    if 1      
        %
        %
        % Reconstruct PASCAL VOC-like synthetic data using various methods
        %
        %
        
        if 1
            % A1. Load pre-rendered synthetic images and just generate some
            % associated metadata
            synth_file = ['synth_5_views_per_mesh'];
            if(~exist([exp_dir '/ImageSets/Segmentation/' synth_file '.txt'], 'file')) % do once
                load_synth_dataset(exp_dir);
            end
        else            
            % A1 (alternative). Generate the synthetic dataset (~3 hours because it first estimates cameras on PASCAL VOC)
            t = tic();
            synth_file = ['synth_' int2str(N_SAMPLES_PER_MESH) '_views_per_mesh'];
            if(~exist([exp_dir '/ImageSets/Segmentation/' synth_file '.txt'], 'file')) % do once
                generate_synth_dataset(exp_dir);
            end
            toc(t)
        end

        % A2. Assemble the full pascal+synth image sets + other associated
        % data    
        if ~exist([exp_dir 'ImageSets/Segmentation/all_gt_segm_kp_' synth_file '.txt'], 'file')
            system(['cat ' exp_dir 'ImageSets/Segmentation/all_gt_segm_kp.txt > ' exp_dir 'ImageSets/Segmentation/all_gt_segm_kp_' synth_file '.txt']);
            system(['cat ' exp_dir 'ImageSets/Segmentation/' synth_file '.txt >> ' exp_dir 'ImageSets/Segmentation/all_gt_segm_kp_' synth_file '.txt']);
            generate_imgset_mirror(exp_dir, ['all_gt_segm_kp_' synth_file]);
            create_mirrored_images(exp_dir, ['all_gt_segm_kp_' synth_file]);
            img_names = textread([exp_dir 'ImageSets/Segmentation/' synth_file '.txt'],'%s');
            for i=1:numel(img_names)
                load([exp_dir 'MySegmentsMat/ground_truth/' img_names{i} '.mat'],'masks');
                masks = {fliplr(masks{1})};
                save([exp_dir 'MySegmentsMat/ground_truth/' img_names{i} '_mirror.mat'],'masks');
            end

            collect_imgset_keypoints(exp_dir,  {['all_gt_segm_kp_' synth_file]});
            collect_imgset_keypoints_mirror(exp_dir,   {['all_gt_segm_kp_' synth_file]}, {['all_gt_segm_kp_' synth_file '_mirror']});
            collect_imgset_keypoints(exp_dir,  {synth_file});
            collect_imgset_keypoints_mirror(exp_dir,  {synth_file}, {[synth_file '_mirror']});
            generate_imgset_mirror(exp_dir, synth_file);

        end
        
        % A3. Reconstruct objects in the synthetic dataset using the PASCAL VOC
        % exemplars. This requires computing synthetic and PASCAL cameras
        % jointly, which takes again ~3 hours.
        eval('config');
        mkdir('Results');
        reconstr_name{1} = ['Reconstructions_Synth_maxdev_' int2str(MAX_DEV) '_nsamples_' int2str(N_SAMPLES)];            

        t = tic();
        n_iter_1 = n_iter_1; % prevent problems with parfor, where used
        n_iter_2 = n_iter_2; 
        N_SAMPLES = N_SAMPLES;
        MAX_DEV = MAX_DEV;
        
        for l=1:20            
            init_rnd_seed();
            reconstruct_GT_getall(l, n_iter_1, n_iter_2, N_SAMPLES, MAX_DEV, reconstr_name{1});
        end
        evaluate_reconstruction_ranking(exp_dir, reconstr_name{1});            
        time_recons_synth = toc(t)

        
        % A4. Same but without the camera refinement stage        
        t = tic();
        reconstr_name{2} = ['Reconstructions_Synth_maxdev_' int2str(MAX_DEV) '_nsamples_' int2str(N_SAMPLES) '_no_refinement'];
        refinement = false;
        parfor l=1:20
            init_rnd_seed();
            reconstruct_GT_getall(l, n_iter_1, n_iter_2, N_SAMPLES, MAX_DEV, reconstr_name{2}, refinement);
        end
        evaluate_reconstruction_ranking(exp_dir, reconstr_name{2});
        time_recons_synth = toc(t)

        
        % A5. With camera refinement stage but without imprinting
        t = tic();
        reconstr_name{3} = ['Reconstructions_Synth_maxdev_' int2str(MAX_DEV) '_nsamples_' int2str(N_SAMPLES) '_no_imprinting'];
        imprinting = false;
        parfor l=1:20
            init_rnd_seed();
            reconstruct_GT_getall(l, n_iter_1, n_iter_2, N_SAMPLES, MAX_DEV, reconstr_name{3}, true, imprinting);
        end
        evaluate_reconstruction_ranking(exp_dir, reconstr_name{3});
        time_recons_synth = toc(t)
        
        
        % A6. Reconstruct using an inflation-based baseline
        t = tic();
        reconstr_name{4} = 'Reconstructions_Synth_Puff_Baseline';
        for l=1:20
            init_rnd_seed();
            reconstruct_GT_puff_baseline(exp_dir, reconstr_name{4}, l);
        end
        evaluate_reconstruction_ranking(exp_dir, reconstr_name{4});
        time_recons_synth_puff = toc(t)
        
        %A7. Baseline reconstruction using the Convex Hull of the rigid shape
        t = tic();
        reconstr_name{5} = 'Reconstructions_Synth_ConvexHull_Baseline';
        parfor l=1:20
            init_rnd_seed();
            reconstruct_GT_convexHull_baseline(l,exp_dir,n_iter_1,n_iter_2,reconstr_name{5});
        end
        evaluate_reconstruction_ranking(exp_dir, reconstr_name{5});
        time_recons_synth_puff = toc(t)

        % just rerun the evaluation functions at the end again 
        % (results are cached)
        fprintf('\n\n\n\n\n');
        for i=1:5
            fprintf('\n\n%s:\n', reconstr_name{i});
            evaluate_reconstruction_ranking(exp_dir, reconstr_name{3});
        end
    end
    
    
    if 1
        % Reconstruct PASCAL VOC
        eval('config');
        mkdir('Results');
        name = 'Pascal';
        
        t = tic();
        parfor l=1:20
            init_rnd_seed();
            reconstruct_pascal_getall(l, n_iter_1, n_iter_2, N_SAMPLES, MAX_DEV, name);
        end
        toc(t)
    end
end


