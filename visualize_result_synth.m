function [] = visualize_result_synth(exp_dir,reconstr_name, sel_class,ind)
%
% example: visualize_result_synth('./VOC/','Reconstructions_Synth_maxdev_15_nsamples_20', 1,1)
%
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
    img_name = sb.img_names{sb.whole_2_img_ids(in_class(reconstructions(1).triples(1)))};
    
    [a,b] = max([reconstructions.score]);    
    r = reconstructions(b);
    

    % load ground truth synthetic mesh
    load([exp_dir 'MyMeshes/ground_truth/' img_name '.mat'], 'tri', 'vertices');
    
    
    im = imread([exp_dir 'JPEGImages/' img_name '.jpg']);

    figure; subplot(1,3,1);
    imshow(im); title('Original Image');
    
    
    subplot(1,3,2);
    trisurf(tri',vertices(1,:)',vertices(2,:),vertices(3,:)); axis equal
    title('Ground truth mesh'); view([0 0 1]); set(gca,'YDir','rev')
    set(gca,'Zdir','rev');
    
    
    subplot(1,3,3);
    trisurf(r.faces,r.vertices(1,:),r.vertices(2,:),r.vertices(3,:)); axis equal
    title('Carvi result');     view([0 0 1]);set(gca,'YDir','rev')
    set(gca,'Zdir','rev');
    
end    
    
    
end

