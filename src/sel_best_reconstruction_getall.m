function [fv2,statistics] = sel_best_reconstruction_getall(all_triples,R,T,kp,really_all_masks, N_VOXELS,axis_masks,flags)

%flags is a structure with the following fields:
%
%flags.is_articulated
%flags.rot_symmetry
%flags.rot_axis (only important if rot_symmetry == 1)
%flags.angle_step (only important if rot_symmetry == 1)

if(sum(abs(R(:,:,all_triples(1,1)))) == 0)
    % If the rotation wasn't estimated.
    fv2.vertices = [];
    fv2.faces = [];
    fv2.triples = [];
    
    return
end

n_images = size(R,3)/2;
if ~flags.is_articulated
    all_triples = [all_triples all_triples+n_images];
end

thisR = R(:,:,all_triples(1,1));
thisT = T(:,all_triples(1,1));
mask = really_all_masks{all_triples(1,1)};
thisKP = kp{all_triples(1,1)};

old_triples = all_triples;

if flags.rot_symmetry
    selR = thisR;
    old_triples = old_triples(1);
    ref_id = 1;
    if(isfield(flags,'angle_step'))
        angle_step = flags.angle_step;
    else
        angle_step = pi/4;
    end
    for th = angle_step:angle_step:(2*pi - angle_step);
        switch flags.rot_axis
            case 'x'
                r1 = [1 0 0;0 cos(th) -sin(th); 0 sin(th) cos(th)];
            case 'y'
                r1 = [cos(th) 0 sin(th);0 1 0; -sin(th) 0 cos(th) ];
            case 'z'
                r1 = [cos(th) -sin(th) 0 ; sin(th) cos(th) 0; 0 0 1];
        end
        r1 = thisR*r1;
        selR = cat(3,selR,r1);
    end
    
    selT = repmat(thisT, [1 size(selR,3)]);
    all_triples = 1:size(selR,3);
    selMasks = cell(size(all_triples));
    for k =1:size(selR,3)
        selMasks{k} = mask;
    end
    [Xvol, Yvol, Zvol, dt,vol] = precompute_data(selR,selT, selMasks , N_VOXELS, ref_id);    
else    
    un_ids = unique(all_triples(:));
    un_ids(un_ids==0) = [];

    ind = zeros(max(un_ids),1);
    ind(un_ids) = 1:length(un_ids);

    all_triples(all_triples~=0) = ind(all_triples(all_triples~=0));

    ref_id = all_triples(1);
    
    [Xvol, Yvol, Zvol, dt,vol] = precompute_data(R(:,:,un_ids), T(:,un_ids), really_all_masks(un_ids), N_VOXELS,ref_id);
end

dt = single(dt);
siz = size(Xvol);

n_triples = size(all_triples,1);

not_projected = ones(n_triples,1)*size(thisKP,1);
dist_proj = 10*ones(n_triples,1);
area_beforeFilling = zeros(n_triples,1);

all_reconstr = repmat(struct('vertices', single(0), 'faces', int32(0)), n_triples, 1);

voxel_bondary = true(siz);
voxel_bondary(2:end-1,2:end-1,2:end-1) = false;

%set to 1 to plot all of the reconstructions.
plot_debug = 0; 

if(plot_debug)
    scrsz = get(0,'ScreenSize');
    h = figure(100);
    set(h,'Position',[1 1 scrsz(3) scrsz(4)]);
end

all_dt_models = single(zeros(size(dt,1),size(all_triples,1)));

for k =1:size(all_triples,1)
    
    ids_sel = all_triples(k,:);
    ids_sel(ids_sel == 0) = [];
     
    volRec = max(dt(:,ids_sel),[],2);    
    
    if(flags.imprinted)
        [volRec,area_beforeFilling(k)] = imprinted_VisualHull(vol,volRec,mask,thisR,thisT,siz);
    end    

    volRec = 1-(2*(volRec<0));    
    volRec = reshape(volRec,siz);
    volRec(voxel_bondary) = abs(volRec(voxel_bondary));
    
    all_dt_models(:,k) = volRec(:);    
    fv2 = isosurface(Xvol,Yvol,Zvol,volRec,0);
    
    if(isempty(fv2.faces))
        disp('empty volume!')
        continue;
    end

    new_vertices = bsxfun(@plus,thisR*(fv2.vertices)',thisT)';

    not_projected(k) = get_points_reprojection(fv2.faces,new_vertices,thisKP);
    dist_proj(k) = main_projections(fv2.faces,fv2.vertices,axis_masks);    

    all_reconstr(k).vertices = single(fv2.vertices);
    all_reconstr(k).faces = int32(fv2.faces);      
    
    if(plot_debug)
        %Plotting stuff
        figure(100);
        subplot(4,5,k);
        hold all;
        
        trisurf(fv2.faces,new_vertices(:,1),new_vertices(:,2),new_vertices(:,3),'FaceColor','red');
        axis equal;
        axis tight;
        view([0 0 1])
        set(gca,'ydir','rev');
        set(gca,'zdir','rev');
        title(num2str(k));
        axis off;              
    end
end

fv2 = all_reconstr;

[scores] = get_ranking(not_projected, dist_proj);

statistics.not_projected = not_projected;
statistics.dist_proj = dist_proj;
statistics.scores = scores;
statistics.area_beforeFilling = area_beforeFilling;


for i=1:numel(fv2)   
    fv2(i).score = scores(i);
    
    fv2(i).dt = all_dt_models(:,i);
    fv2(i).triples = old_triples(i,:);
end
        
end


function [volRec,area_beforeFilling] = imprinted_VisualHull(vol,volRec,mask,thisR,thisT,siz)

vol = bsxfun(@plus,thisR*vol,thisT);

[M,N] = size(mask);

x = vol(1,:);
y = vol(2,:);

ma = false(size(x));

sel_ind = x >0.5 & x<N & y>0.5 & y<M;
ind = sub2ind(size(mask),round(y(sel_ind)),round(x(sel_ind))); 
ma(sel_ind) = mask(ind);

volRec = reshape(volRec,siz);

min_voxel = repmat(min(volRec,[],3),[1 1 size(volRec,3)]);

ma = ma(:)';
ma = reshape(ma,siz);

area_beforeFilling = sum(sum(ma(:,:,1) &  min_voxel(:,:,1)<0))/sum(sum(ma(:,:,1)));

volRec(ma & min_voxel>0 ) = volRec(ma & min_voxel>0) - (min_voxel(ma & min_voxel>0) + 0.0001);
volRec = volRec(:);

end


function dist_main = main_projections(tri,coord,axis_masks)
% returns the distance between the projections into the planes aligned with the main axis and the
% average mask for those planes

m = axis_masks.size;

edge = zeros(2,0,'uint32');

coord = coord*axis_masks.scale;
coord = coord + m/2;

projections = zeros(m,m,3);


P = [ 0 0 1 0; 0 1 0 0 ; 0 0 0 1];
projections(:,:,1) = RenderTriMex(P, m, m, coord', edge, uint32(tri-1)')>0;

P = [ 0 0 1 0;1 0 0 0 ; 0 0 0 1];
projections(:,:,2) = RenderTriMex(P, m, m, coord', edge, uint32(tri-1)')>0;

P = [ 0 1 0 0; 1 0 0 0; 0 0 0 1];
projections(:,:,3) = RenderTriMex(P, m, m, coord', edge, uint32(tri-1)')>0;

dist_main = 0;

for k =1:3
    if(axis_masks.exemplars(k)~=0);
        dist_main = dist_main + sum(sum(abs(projections(:,:,k) - axis_masks.masks(:,:,k))));
    end
end

dist_main = dist_main/(m*m);


end


function [not_projected,thisKP_3D] = get_points_reprojection(tri,rot_vertices,thisKP)
%Returns the number of points that does not reproject inside the shape

TR2D = TriRep(tri,rot_vertices(:,1:2));
TR3D = TriRep(tri,rot_vertices);

n_tri = size(tri,1);

not_projected = 0;
thisKP_3D = zeros(size(thisKP,1),3);
for p = 1:size(thisKP,1)
    if(thisKP(p,1)==0 && thisKP(p,2)==0); continue; end;

    B = cartToBary(TR2D,(1:n_tri)',repmat(thisKP(p,:),[n_tri 1]));

    sel_triangles = B>=0 & B<=1;

    sel_triangles = find(sum(sel_triangles,2) == 3);
    if(~isempty(sel_triangles))            
        PC = baryToCart(TR3D,sel_triangles,B(sel_triangles,:));
  
        [a,b] = min(PC(:,3));
        thisKP_3D(p,:) = PC(b,:);
    else
        not_projected = not_projected+1;
    end

end
        
end



function [Xvol, Yvol, Zvol, dt, vol] = precompute_data(R, T, really_all_masks, n_voxels,ref_id)
    
    lim = get_vol_limits(really_all_masks, R, T,ref_id);
    
    for i=1:3
        step(i) = (lim(i,2) - lim(i,1)) / n_voxels;
    end
    
    [Xvol,Yvol,Zvol] = meshgrid(lim(1,1):step(1):lim(1,2),lim(2,1):step(2):lim(2,2),lim(3,1):step(3):lim(3,2));    
    sz = size(Xvol);
    vol = [Xvol(:)'; Yvol(:)'; Zvol(:)'];
    
    refR = R(:,:,ref_id);
    refT = T(:,ref_id);
    
    vol = refR\bsxfun(@minus,vol, refT);
    
    Xvol = reshape(vol(1,:)',sz);
    Yvol = reshape(vol(2,:)',sz);
    Zvol = reshape(vol(3,:)',sz);

    
    dt = zeros(size(vol,2),size(R,3));
    
    for h=1:size(R,3)
        thisR = R(:,:,h);
        thisT = T(:,h);
        mask = really_all_masks{h};
        
        rotVol = thisR*vol;
        rotVol =  bsxfun(@plus,rotVol,thisT);
        dt(:,h) = distance_from_mask(mask,rotVol);
    end
end


function dt = distance_from_mask(mask,rotVol)
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


dt = zeros(size(rotVol,2),1);

[x,y] = getXYfromMask(aux_mask);

[~,IDX] = bwdist(aux_mask);
sel_points = IDX(ind);
d = sqrt((x(sel_points)-rotVol(1,:)').^2 + (y(sel_points)-rotVol(2,:)').^2);

dt(~aux_mask(ind)) = d(~aux_mask(ind));

[~,IDX] = bwdist(~aux_mask);
sel_points = IDX(ind);
d = sqrt((x(sel_points)-rotVol(1,:)').^2 + (y(sel_points)-rotVol(2,:)').^2);

dt(aux_mask(ind)) = -d(aux_mask(ind));

end


function [x,y,m] = getXYfromMask(mask)

[x,y] = meshgrid(1:size(mask,2),1:size(mask,1));

x = x(:);
y = y(:);
m = mask(:);


end



function [limits] = get_vol_limits(masks, R, T,ref_id)   


refR = R(:,:,ref_id);
refT = T(:,ref_id);


min_x = inf;
min_y = min_x;
min_z = min_x;
max_x = -inf;
max_y = max_x;
max_z = max_x;

for h=1:size(R,3)
    % get edgels
    bw = bwboundaries(masks{h});
    bw = bw{1};


    x = bw(:,2);
    y = bw(:,1);
    z = zeros(size(bw,1),1);

    coord = [x'; y';z'];

    thisR = R(:,:,h);
    thisT = T(:,h);


    coord = thisR\(bsxfun(@minus, coord, thisT));
    coord = bsxfun(@plus, refR*coord, refT);

    min_x = min([min_x coord(1,:)]);
    max_x = max([max_x coord(1,:)]);

    min_y = min([min_y coord(2,:)]);
    max_y = max([max_y coord(2,:)]);

    min_z = min([min_z coord(3,:)]);
    max_z = max([max_z coord(3,:)]);
end

limits = [min_x max_x; min_y max_y; min_z max_z];
end




