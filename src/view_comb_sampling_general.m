function [ids_target, ids_pool, ax] = view_comb_sampling_general(R_target, R_pool, N_SAMPLES_PER_OBJ, max_angle_deviation)
   % 
   DefaultVal('*N_SAMPLES_PER_OBJ', '20');
   DefaultVal('*max_angle_deviation', '[40 20 25]'); % in degrees
   
    % for each image, sample another two images belonging to one of the
    % three clusters of images orthogonal to the objects left-right
    % symmetry axis.

    % 1. identify sets of axis aligned images 
    n_obj_target = size(R_target,3)/2;
    [ax_target{1}, ax_target{2}, ax_target{3}] = collect_axis(R_target(:,:,1:n_obj_target), max_angle_deviation);

    ax{1} = [];
    ax{2} = [];
    ax{3} = [];
    % only one axis is allowed to be empty and the others should have at
    % least 2 elements so that if our test view is one of them we don't get
    % a second empty axis (sofa was crashing)
    n_ex = sort(cellfun(@numel, ax),'ascend');
    while(n_ex(2)<2)
        n_obj_pool = ceil(size(R_pool,3)/2);
        [ax{1}, ax{2}, ax{3}] = collect_axis(R_pool(:,:,1:n_obj_pool), max_angle_deviation);

        % each object should belong to a single axis, if any
        assert(isempty(intersect(ax{1},ax{2})));
        assert(isempty(intersect(ax{2},ax{3})));
        assert(isempty(intersect(ax{1},ax{3})));

        max_angle_deviation = max_angle_deviation + 5;
        n_ex = sort(cellfun(@numel, ax),'ascend');
    end

    if 0
        if(~isempty(ax{1}))
            figure; sb.show_wholes(in_class(ax{1}));
        end
        if(~isempty(ax{2}))
            figure;sb.show_wholes(in_class(ax{2}));
        end
        if(~isempty(ax{3}))
            figure;sb.show_wholes(in_class(ax{3}));
        end
    end

    ids_target = [];
    ids_pool = [];
    
    n_exemplars = cellfun(@(x)length(x),ax);


    for i=1:n_obj_target
        in_axis = [];
        if any(ax_target{1} == i)
            in_axis = 1;
        elseif any(ax_target{2} == i)
            in_axis = 2;
        elseif any(ax_target{3} == i)
            in_axis = 3;
        end

        %probability of selecting each of the axis
        % It should take into account:
        % - number of examples in that axis
        % - if the object belongs to that axis
        % - it should be zero if there is no element
        axis_probability = n_exemplars./sum(n_exemplars);
        axis_probability = max(axis_probability,0.2);
        if(n_exemplars(in_axis) == 1)
            axis_probability(in_axis) = 0;
        else
            axis_probability(in_axis) = 0.1;
        end
        axis_probability(n_exemplars == 0) = 0;
        a = 1-sum(axis_probability(axis_probability==0.1));
        n = sum(axis_probability(axis_probability ~= 0.1));
        axis_probability(axis_probability~=0.1) = axis_probability(axis_probability~=0.1)*a/n;
        cum_probability = cumsum(axis_probability);

        rand_ax = rand(N_SAMPLES_PER_OBJ,1);
        sel_axis1 = zeros(N_SAMPLES_PER_OBJ,1);

        for k =1:3
            sel_axis1(rand_ax<=cum_probability(k) & sel_axis1==0) = k;
        end

        axis_probability = repmat(axis_probability,[N_SAMPLES_PER_OBJ 1]);

        ind = sub2ind(size(axis_probability),(1:N_SAMPLES_PER_OBJ)',sel_axis1);
        axis_probability(ind) = 0;

        axis_probability = bsxfun(@rdivide,axis_probability,sum(axis_probability,2));
        cum_probability = cumsum(axis_probability,2);

        rand_ax = rand(N_SAMPLES_PER_OBJ,1);

        sel_axis2 = zeros(N_SAMPLES_PER_OBJ,1);

        for k =1:3
            sel_axis2(rand_ax<=cum_probability(:,k) & sel_axis2==0) = k;
        end

        sel_axis = [sel_axis1 sel_axis2];

        n_sel_ax = hist(sel_axis(:),[1 2 3]);

        second_third = zeros(N_SAMPLES_PER_OBJ,2);

        for k = 1:3
            if in_axis == k
                sel_ids = ceil(rand(n_sel_ax(k),1)*(n_exemplars(k)-1));
                cand = ax{k};
                cand(cand == i) = [];
            else
                sel_ids = ceil(rand(n_sel_ax(k),1)*n_exemplars(k));
                cand = ax{k};
            end
            cand = cand(sel_ids);
            second_third(sel_axis==k) = cand;
        end

        ids_target = [ids_target; i*ones(N_SAMPLES_PER_OBJ,1)];
        ids_pool = [ids_pool; second_third];        

        % Just the image itself. When alignment is hard the image and its symmetric may lead
        % to better reconstruction than pairs/triplets
        % but it has to be carefully taken when evaluating...
        %ids = [ids; i 0 0];
    end
end
