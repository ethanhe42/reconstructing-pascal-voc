function [ axis_masks ] = get_average_mask(axs, R, T, really_all_masks)

axis_inds = [1 2 3];

axis_masks.masks = zeros(200,200,3);

axis_masks.scale = 3;
axis_masks.size = 200;
axis_masks.exemplars = [ 0 0 0];

for i = axis_inds

    sel = [1 2 3];
    sel(sel ==i) = [];
    ax = axs{i};

    sum_mask = zeros(200);

    for k = ax(1:end)
        thisR = R(:,:,k);
        sc = norm(thisR(1,:))/4;
        
        if(sc<0.1)
            continue;
        end
        
        thisT = T(:,k);
        mask = really_all_masks{k};

        [M,N] = size(mask);
        
        [x,y] = meshgrid(1:sc:N,1:sc:M);
        mask = imresize(mask,size(x));
        x = x(:)';
        y = y(:)';

        coord = [x;y;zeros(size(x))];
        coord = coord(:,mask(:));

        coord = bsxfun(@minus,coord,thisT);
        coord = axis_masks.scale*(thisR\coord);
        coord = coord(sel,:);

        coord = bsxfun(@plus,coord,[axis_masks.size/2; axis_masks.size/2]);

        coord(1,:) = max(1,min(coord(1,:),size(sum_mask,2)));
        coord(2,:) = max(1,min(coord(2,:),size(sum_mask,1)));


        ind = sub2ind(size(sum_mask),round(coord(2,:)'),round(coord(1,:)'));
        B = zeros(size(sum_mask));
        B(ind) = 1;
        sum_mask(ind) = sum_mask(ind) +1;

    end

    if(length(ax)~=0)
        sum_mask = sum_mask/length(ax);
    end
    
    axis_masks.exemplars(i) = length(ax);
    axis_masks.masks(:,:,i) = sum_mask; 
end




end

