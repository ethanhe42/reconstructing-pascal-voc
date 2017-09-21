function [axis_1, axis_2, axis_3] = collect_axis(R, max_dev_angle, view_based)
    DefaultVal('*view_based', 'false');
    p = [0 0 1]';
    a1 = [1 0 0]';
    a2 = [0 1 0]';
    a3 = p;

    n_images = size(R,3);
    dir_vector = zeros(3,n_images);
    for k =1:n_images
        thisR = R(:,:,k);
        s = norm(thisR(1,:));
        if(s <0.01)
            dir_vector(:,k) = [nan nan nan];
        else
            thisR = thisR/s;
            invR = thisR';
            dir_vector(:,k) = invR*p;
        end
        if(any(isnan(dir_vector(:,k))))
            %disp('bad rotation matrix');
            %R(:,:,k)
        end
    end

    
    a = bsxfun(@times,dir_vector,a1);
    cos_angles_1 = sum(a);
        
    a = bsxfun(@times,dir_vector,a2);
    cos_angles_2 = sum(a);
        
    a = bsxfun(@times,dir_vector,a3);
    cos_angles_3 = sum(a);
    
    %cos_angles_1(isnan(cos_angles)) = inf; % set high so we can later reject it
    
    the_axis{1} = find(abs(cos_angles_1) > cos(deg2rad(max_dev_angle(1)))); % radians
    the_axis{2} = find(abs(cos_angles_2) > cos(deg2rad(max_dev_angle(2))));
    the_axis{3} = find(abs(cos_angles_3) > cos(deg2rad(max_dev_angle(3))));
    
    if 0 % (view_based)
        % copy the axis with most elements to the axis with fewest (view
        % based means there is no side view, so that one will be empty and
        % the front view will be abundant)

        n = cellfun(@numel, the_axis);
        [val_max,id_max] = max(n);
        [val_min,id_min] = min(n);
        
        the_axis{id_min} = [the_axis{id_min} the_axis{id_max}];
                
        % gotta_go = setdiff([1 2 3], [id_min id_max]);
        % the_axis{gotta_go} = [];
    end
    
    axis_1 = the_axis{1};
    axis_2 = the_axis{2};
    axis_3 = the_axis{3};
end