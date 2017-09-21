function reconstruct_GT_puff_baseline_pascal(sel_class)
    OFFICE = true;
    exp_dir = add_all_paths(OFFICE);
    
    DefaultVal('*sel_class', '1');    

    mask_type = 'ground_truth';       

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% Get data of images having ground truth meshes %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    imgset_all = 'all_gt_segm_kp';
    pascal_mask = 'ground_truth';
    sb_pascal = SegmBrowser(exp_dir, pascal_mask, imgset_all);
    y = sb_pascal.get_overlaps_wholes(1:numel(sb_pascal.whole_2_img_ids));
    [q, class] = max(y, [], 2);
    in_class = find(class == sel_class);    
    masks = sb_pascal.get_masks(in_class);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Reconstruct object instances %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %all_reconstr = repmat(struct('vertices', single(0), 'faces', int32(0),'dt', single(0)), size(in_class,1), 1);
    t_total = tic();
    
    out_folder = ['SFM_models/Reconstructions_pascal_Puff_Baseline/' sb_pascal.categories{sel_class} '/'];
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
            I = sb_pascal.get_Imgs(sb_pascal.whole_2_img_ids(in_class));
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

function [Xvol, Yvol, Zvol, dt] = precompute_data(R, T, Shape, really_all_masks, step, border, faster)
    %TODO: this should probably be adjusted depending on the class, to ensure
    %that the volumetric representation covers all the volume it needs.    
    mi = floor(min(Shape,[],2)) -border;
    ma = ceil(max(Shape,[],2)) + border;        
    [Xvol,Yvol,Zvol] = meshgrid(mi(1):step:ma(1),mi(2):step:ma(2),mi(3):step:ma(3));    
    vol = [Xvol(:)'; Yvol(:)'; Zvol(:)'];
    dt = zeros(size(vol,2),size(R,3));
    
    %parfor h=1:size(R,3)
    for h=1:size(R,3)
        thisR = R(:,:,h);
        thisT = T(:,h);
        mask = really_all_masks{h};
        
        rotVol = thisR*vol;
        rotVol =  bsxfun(@plus,rotVol,thisT);
        dt(:,h) = distance_from_mask(mask,rotVol,faster);
    end
end

function dt = distance_from_mask(mask,rotVol,faster)
% Returns a negative value for voxels in the volume and positive value
% outside.


[M,N] = size(mask);

minX = floor(min(rotVol(1,:))-2);
minY = floor(min(rotVol(2,:))-2);

maxX = ceil(max(rotVol(1,:))+2);
maxY = ceil(max(rotVol(2,:))+2);

cX = max(0,-minX+1);
cY = max(0,-minY+1);

aux_mask = false(cY + max(maxY,M),cX + max(maxX,N));

aux_mask(cY+1:cY+M, cX+1:cX+N) = mask;

rotVol(1,:) = rotVol(1,:) + cX;
rotVol(2,:) = rotVol(2,:) + cY;

ind = sub2ind(size(aux_mask),round(rotVol(2,:)),round(rotVol(1,:)));

if(faster)
    dt = -2*double(aux_mask(ind))' + 1;
else
    dt = zeros(size(rotVol,2),1);

    [x,y] = getXYfromMask(aux_mask);
    
    %Discrete approximation of distance, it simply takes the closest point in
    %terms of the 2D distance transform, and then refines that distance.
    [~,IDX] = bwdist(aux_mask);
    sel_points = IDX(ind);
    d = sqrt((x(sel_points)-rotVol(1,:)').^2 + (y(sel_points)-rotVol(2,:)').^2);
    
    dt(~aux_mask(ind)) = d(~aux_mask(ind));
    
    [~,IDX] = bwdist(~aux_mask);
    sel_points = IDX(ind);
    d = sqrt((x(sel_points)-rotVol(1,:)').^2 + (y(sel_points)-rotVol(2,:)').^2);
    
    dt(aux_mask(ind)) = -d(aux_mask(ind));
end

end


function [x,y,m] = getXYfromMask(mask)

[x,y] = meshgrid(1:size(mask,2),1:size(mask,1));

x = x(:);
y = y(:);
m = mask(:);


end

function out = assig2cell(assig)
    un = unique(assig);
    n = numel(un);
    
    out = cell(n,1);
    for i=1:n
        out{i} = find(assig==un(i));
    end
    out(isempty(out)) = [];
end