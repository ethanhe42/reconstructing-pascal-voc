function load_synth_dataset(exp_dir)
    load('voc_kp_metadata.mat');
    eval('config');
    
    %%%%%%%%%%%%%%%%%%% params & directories & filenames %%%%%%%%%%%%%%%%%%            
    imgset_name = 'synth_5_views_per_mesh';    
     
    system(['cp ./Dataset/Renderings/JPEGImages/* ' exp_dir 'JPEGImages/']);    
    
    system(['cp ./Dataset/Renderings/MySegmentsMat/ground_truth/* ' exp_dir 'MySegmentsMat/ground_truth/']);
            
    mkdir([exp_dir 'merged_Correspondences_GT_BRKL/']);
    system(['cp ./Dataset/Renderings/merged_Correspondences_GT_BRKL/* ' exp_dir 'merged_Correspondences_GT_BRKL/']);
                    
    system(['cp ./Dataset/Renderings/' imgset_name '.txt ' exp_dir '/ImageSets/Segmentation/']);    
            
    system(['cp ./Dataset/Renderings/merged_Correspondences_GT_BRKL/* ' exp_dir 'merged_Correspondences_GT_BRKL/']);

    system(['cp ./Dataset/Renderings/SegmentationObject/* ' exp_dir '/SegmentationObject/']);
    system(['cp ./Dataset/Renderings/SegmentationClass/* ' exp_dir '/SegmentationClass/']);
    
    mkdir([exp_dir 'MyMeshes/ground_truth/']);
    system(['cp ./Dataset/Renderings/MyMeshes/ground_truth/* ' exp_dir 'MyMeshes/ground_truth/']);
        
    for i=1:20
        categs{i} = VOC09_id_to_classname(i);
    end
    
    for h=1:numel(categs)        
        mesh_files = dir(['./Dataset/SynthMeshes/' categs{h} '/*_mesh.mat']);
        n_meshes = numel(mesh_files);
                
        for i=1:n_meshes            
            % 5 examples of each of 10 meshes per class
            base_output_img_name = strtok(mesh_files(i).name,'_');

            for j=1:5
                output_img_name = [base_output_img_name '_view_' int2str(j)];
                
                var = load(['./Dataset/SynthMeshes/' categs{h} '/' base_output_img_name '_mesh.mat']);
                tri = var.nodes;
                
                load([exp_dir 'MyMeshes/ground_truth/' output_img_name]);

                transf_xyz = bsxfun(@plus, global_s*global_R*var.xyz, global_T);
                cam_xyz = bsxfun(@plus,camR*transf_xyz,camT);
                vertices = cam_xyz;
                
                save([exp_dir 'MyMeshes/ground_truth/' output_img_name], 'tri', 'vertices', 'global_R','global_s', 'global_T', 'P', 'cam_transf_kp', 'camR', 'camT', 'visible_kp');
            end
        end
    end
    
end