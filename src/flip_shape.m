function [R,Shape] = flip_shape(R,Shape,all_kp,right_coord_sys,sel_kp_ids)

n_im = size(R,3);

vis = zeros(n_im,1);
not_vis = zeros(n_im,1);


for k = 1:n_im
    kp = all_kp(:,:,k);
    not_visible = isnan(kp(1,:));
    
    rot_shape = R(:,:,k)*Shape;
    
    vis(k) = sum(rot_shape(3,~not_visible));
    
    not_vis(k) = sum(rot_shape(3,not_visible));
    
end

th = sum(vis>not_vis)/numel(vis);
if(abs(th-0.5) <0.2 && ~isempty(right_coord_sys))
    aux_ids = zeros(max(max(sel_kp_ids(:)),max(right_coord_sys(:))),1);
    aux_ids(sel_kp_ids) = 1:numel(sel_kp_ids);
    
    right_coord_sys = aux_ids(right_coord_sys);
    
    a = Shape(:,right_coord_sys(1,2)) - Shape(:,right_coord_sys(1,1));
    b = Shape(:,right_coord_sys(2,2)) - Shape(:,right_coord_sys(2,1));
    c = cross(a,b);
    v = Shape(:,right_coord_sys(3,2)) - Shape(:,right_coord_sys(3,1));
    
    if(dot(c,v)<0)
        disp('flipped');

        Shape = -Shape;

        R(1:2,:,:) = -R(1:2,:,:);
    end
elseif(th>0.5)
    disp('flipped');
    
    Shape = -Shape;
    
    R(1:2,:,:) = -R(1:2,:,:);
end


end
