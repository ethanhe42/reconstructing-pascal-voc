function assemble_voc_kp_metadata()    
    exp_dir = add_all_paths();
    
    mask_type_train = 'ground_truth';

    imgset_all = 'all_gt_segm';

    sb = SegmBrowser(exp_dir, mask_type_train, imgset_all);
    corr = CorrespBRKL(exp_dir, imgset_all, sb);
    
    metadata.categories = sb.categories;
    metadata.articulated = logical([0 0 1 1 0 0 0 1 0 1 0 1 1 0 1 0 1 0 0 0]); % bicycle and train approximated as rigid, which is true most of the time 
    metadata.view_based = logical([0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0]);
    
    for i=1:numel(sb.categories)
        %%%%%%%%%%%% Keypoint names %%%%%%%%%%%%%%
        metadata.kp_names{i} = corr.get_category_keypoint_names(i);
                
        %%%%%%%%%%%% Rigidity   %%%%%%%%%%%%%%%%%%
        if 1
            % these are used for selecting keypoints to do rigid SFM over
            % (by selecting largest rigid components)
            metadata.rigidity{i} = false(numel(metadata.kp_names{i}));
            if(~metadata.articulated(i))
                metadata.rigidity{i} = true(numel(metadata.kp_names{i})); % everything is rigidly connected
                metadata.rigid_parts{i}{1} = 1:numel(metadata.kp_names{i});
            else
                switch sb.categories{i}                    
                    case 'bird'                      
                        ids = [1:3 8];
                        metadata.rigid_parts{i}{1} = ids;
                        metadata.rigidity{i}(ids, ids) = true; % head
                        
                        ids = [4 5];
                        metadata.rigid_parts{i}{2} = ids;
                        metadata.rigidity{i}(ids, ids) = true; % left wing
                        
                        ids = [9 10];
                        metadata.rigid_parts{i}{3} = ids;
                        metadata.rigidity{i}(ids, ids) = true; % right wing
                                                
                        ids = [4 9 6 7 11 12];
                        metadata.rigid_parts{i}{4} = ids;
                        metadata.rigidity{i}(ids, ids) = true; % torso    
                    case 'boat'
                        ids = [1:8];
                        metadata.rigid_parts{i}{1} = ids;
                        metadata.rigidity{i}(ids, ids) = true; % hull
                        
                        ids = [9:11];
                        metadata.rigid_parts{i}{2} = ids;
                        metadata.rigidity{i}(ids, ids) = true; % sail                        
                    case 'cat'                       
                        ids = [3 4 7 10 11 15];
                        metadata.rigid_parts{i}{1} = ids;
                        metadata.rigidity{i}(ids, ids) = true; % head
                                                
                        ids = [12 13];
                        metadata.rigid_parts{i}{2} = ids;
                        metadata.rigidity{i}(ids, ids) = true; % right frontal leg
                        
                        ids = [8 9];
                        metadata.rigid_parts{i}{3} = ids;
                        metadata.rigidity{i}(ids, ids) = true; % right rear leg
                        
                        ids = [5 6];
                        metadata.rigid_parts{i}{4} = ids;
                        metadata.rigidity{i}(ids, ids) = true; % left frontal leg
                        
                        ids = [1 2];
                        metadata.rigid_parts{i}{5} = ids;
                        metadata.rigidity{i}(ids, ids) = true; % left rear leg
                        
                        ids = [1 5 8 12 14 15 16];
                        metadata.rigid_parts{i}{6} = ids;
                        metadata.rigidity{i}(ids, ids) = true; % torso                                                
                    case 'cow'
                        metadata.rigidity{i} = metadata.rigidity{8}; % same as cat                        
                        metadata.rigid_parts{i} = metadata.rigid_parts{8};
                    case 'dog'
                        metadata.rigidity{i} = metadata.rigidity{8}; % same as cat  
                        metadata.rigid_parts{i} = metadata.rigid_parts{8};
                    case 'horse'
                        metadata.rigidity{i} = metadata.rigidity{8}; % same as cat  
                        metadata.rigid_parts{i} = metadata.rigid_parts{8};
                    case 'person'                                              
                       ids = [1 2 4 6 13 15 17];
                       metadata.rigid_parts{i}{1} = ids;
                       metadata.rigidity{i}(ids, ids) = true; % head
                       
                       ids = [8 19 10 21];
                       metadata.rigid_parts{i}{2} = ids;
                       metadata.rigidity{i}(ids, ids) = true; % torso
                       
                       ids = [3 7 9 11];
                       metadata.rigid_parts{i}{3} = ids;
                       metadata.rigidity{i}(ids, ids) = true; % lower left leg
                       
                       ids = [8 9];
                       metadata.rigid_parts{i}{4} = ids;
                       metadata.rigidity{i}(ids, ids) = true; % upper left leg
                       
                       ids = [14 18 20 22];
                       metadata.rigid_parts{i}{5} = ids;
                       metadata.rigidity{i}(ids, ids) = true; % lower right leg
                       
                       ids = [19 20];
                       metadata.rigid_parts{i}{6} = ids;
                       metadata.rigidity{i}(ids, ids) = true; % upper right leg                       
                       
                       ids = [5 12];
                       metadata.rigid_parts{i}{7} = ids;
                       metadata.rigidity{i}(ids, ids) = true; % lower left arm
                       
                       ids = [10 5];
                       metadata.rigid_parts{i}{8} = ids;
                       metadata.rigidity{i}(ids, ids) = true; % upper left arm
                       
                       ids = [16 23];
                       metadata.rigid_parts{i}{9} = ids;
                       metadata.rigidity{i}(ids, ids) = true; % lower right arm
                       
                       ids = [21 16];
                       metadata.rigid_parts{i}{10} = ids;
                       metadata.rigidity{i}(ids, ids) = true; % upper right arm                      
                    case 'sheep'
                        metadata.rigidity{i} = metadata.rigidity{8}; % same as cat             
                        metadata.rigid_parts{i} = metadata.rigid_parts{8};
                    otherwise
                        error('no such other type');                
                end
            end
        end
        
        %%%%%%%%%%%% Left-Right Symmetries %%%%%%%%%%%%%%%%%%
        if 1
            %corr.show_categ_objects(i, [1:5]);

            switch sb.categories{i}                
                case 'aeroplane'
                    metadata.sym_corresp{i} = {[3 4 5 6 7; 11 12 13 14 15]};
                case 'bicycle'
                    metadata.sym_corresp{i} = {[9; 10]};
                case 'bird'
                    metadata.sym_corresp{i} = {[3 4 5; 8 9 10]};
                case 'boat'
                    metadata.sym_corresp{i} = {[5 6 10; 7 8 11]}; 
                case 'bottle'
                    metadata.sym_corresp{i} = {[1:4; 5:8]}; 
                case 'bus'
                    metadata.sym_corresp{i} = {[1:4; 5:8]}; 
                case 'car'
                    metadata.sym_corresp{i} = {[1:7; 8:14]};
                case 'cat' 
                    metadata.sym_corresp{i} = {[1 2 3 4 5 6; 8 9 10 11 12 13]}; 
                case 'chair' 
                    metadata.sym_corresp{i} = {[1 3 4 7 8; 2 5 6 9 10]};
                case 'cow' 
                    metadata.sym_corresp{i} = metadata.sym_corresp{8}; % same as cat
                case 'diningtable'
                    metadata.sym_corresp{i} = {[1 2 5 6; 3 4 7 8]};
                case 'dog'
                    metadata.sym_corresp{i} = metadata.sym_corresp{8}; % same as cat
                case 'horse' 
                    metadata.sym_corresp{i} = metadata.sym_corresp{8}; % same as cat
                case 'motorbike'
                    metadata.sym_corresp{i} = {[7; 8]}; % similar to bicycle
                case 'person'
                    metadata.sym_corresp{i} = {[3:12; 14:23]}; 
                case 'pottedplant'
                    metadata.sym_corresp{i} = {[1 5; 2 6]}; 
                case 'sheep'
                    metadata.sym_corresp{i} = metadata.sym_corresp{8}; % same as cat      
                case 'sofa'
                    metadata.sym_corresp{i} = {[1 3 5 7 9 11; 2 4 6 8 10 12]}; 
                case 'train'
                    metadata.sym_corresp{i} = {[1 3 5; 2 4 6]}; 
                case 'tvmonitor'
                    metadata.sym_corresp{i} = {[1 3 5 7; 2 4 6 8]}; 
                otherwise
                    error('no such other type');            
            end
        end  
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% to determin rotation axis for view-based classes %%%%%%%%%%%%%%%
    for i=find(metadata.view_based)
        switch sb.categories{i}
            case 'bottle'
                metadata.bottom_top{i} = {[1 5],[4 8]};
            case 'diningtable'
                metadata.bottom_top{i} = {1:4,5:8};
            case 'pottedplant'
                metadata.bottom_top{i} = {1:2,3:6};
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% To correct flip %%%%%%%%%%%%%%%
    for i=1:numel(sb.categories)
        switch sb.categories{i}
            case 'aeroplane'
                metadata.right_coordinate_sys{i} = [16 11;16 3;16 8];
            case 'bicycle'
                metadata.right_coordinate_sys{i} = [5 10;5 9; 5 1];
            case 'bird'
                metadata.right_coordinate_sys{i} =  [9 12;9 4; 9 6];
            case 'boat'
                metadata.right_coordinate_sys{i} = [6 8; 6 7; 6 3];
            case 'bottle'
                metadata.right_coordinate_sys{i} = [];
            case 'bus'
                metadata.right_coordinate_sys{i} = [4 3; 4 7; 4 2];
            case 'car'
                metadata.right_coordinate_sys{i} = [3 10; 3 6; 3 5];
            case 'cat'
                metadata.right_coordinate_sys{i} = [5 16; 5 12; 5 7];
            case 'chair'
                metadata.right_coordinate_sys{i} =  [3 1;3 5; 3 4];
            case 'cow'
                metadata.right_coordinate_sys{i} = [5 16; 5 12; 5 7];% same as cat
            case 'diningtable'
                metadata.right_coordinate_sys{i} =  [6 8; 6 2; 6 1];
            case 'dog'
                metadata.right_coordinate_sys{i} = [5 16; 5 12; 5 7];% same as cat
            case 'horse'
                metadata.right_coordinate_sys{i} = [5 16; 5 12; 5 7];% same as cat
            case 'motorbike'
                metadata.right_coordinate_sys{i} = [7 4;7 8;7 10];
            case 'person'
                metadata.right_coordinate_sys{i} = [13 21;13 10;13 19];
            case 'pottedplant'
                metadata.right_coordinate_sys{i} =  [5 1;5 6;5 4];
            case 'sheep'
                metadata.right_coordinate_sys{i} = [5 16; 5 12; 5 7];% same as cat
            case 'sofa'
                metadata.right_coordinate_sys{i} =  [2 4; 2 6; 2 1];
            case 'train'
                metadata.right_coordinate_sys{i} = [4 7; 4 3; 4 2];
            case 'tvmonitor'
                metadata.right_coordinate_sys{i} = [7 5;7 8;7 3];
            otherwise
                error('no such other type');
        end
    end
    
    save('voc_kp_metadata.mat', 'metadata');        
end
