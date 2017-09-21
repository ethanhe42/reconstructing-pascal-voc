function [] = visualize_result(exp_dir,reconstr_name, sel_class,ind)
    %
    % example: visualize_result('./VOC/','Pascal', 9,6)
    %

    add_all_paths();
    load voc_kp_metadata

    in_class = [];

    base_dir = ['./Results/' reconstr_name '/' metadata.categories{sel_class} '/'];
    load([base_dir 'reconstruction_data.mat']);

    sb = SegmBrowser(exp_dir, 'ground_truth', imgset_all);

    files = dir([base_dir '*.mat']);
    files(arrayfun(@(a) strcmp(a.name, 'reconstruction_data.mat'), files)) = [];

    if(nargin ~= 4 || isempty(ind))
        ind = randperm(length(files),1);
    end

    for k = ind
        disp(k);
        % load reconstructed mesh
        load([base_dir files(k).name], 'reconstructions');

        [a,b] = max([reconstructions.score]);    
        r = reconstructions(b);       

        figure;
        im = sb.show_wholes(in_class(reconstructions(1).triples(1)));
        subplot(1,2,1);
        imshow(im{1}); title('Original Image');

        subplot(1,2,2);
        trisurf(r.faces,r.vertices(1,:),r.vertices(2,:),r.vertices(3,:)); axis equal
        title('Carvi result');     view([0 0 1]);set(gca,'YDir','rev')
        set(gca,'Zdir','rev');    
    end    
end

