%function SvmSegm_study_segment_quality(exp_dir, segm_name, the_imgset, segm_quality_type)
function SvmSegm_study_segment_quality(exp_dir, segm_name, the_imgset, segm_quality_type, class_label)
    DefaultVal('class_label', '[]');
       
    if(~iscell(segm_name))
        segm_name = {segm_name};
    end

    % can write segment statistics into [exp MyStatistics]
    % one file per segmentation algorithm
    %statistics_dir = [exp_dir 'MyStatistics/'];
    %if(~exist(statistics_dir, 'dir'))
    %    mkdir(statistics_dir);
    %end

    if(~iscell(the_imgset))
        img_names =textread([exp_dir '/ImageSets/Segmentation/' the_imgset '.txt'], '%s');
    else
        if(5 > numel(the_imgset))
            error('if it''s a imgset, don''t pass as cell!');
        end
        img_names = the_imgset;
    end

    segm_eval_dir = [exp_dir 'SegmentEval/'];

    Q = {};
            
    for i=1:numel(segm_name)
        mask_type = segm_name;
        direct = [segm_eval_dir mask_type{i} '/' segm_quality_type '/' ];
        if(~exist(direct, 'dir'))
            error('directory doesn''t exist');
        end
        
        thisQ = cell(numel(img_names),1);
        for j=1:numel(img_names)
            if(exist([direct img_names{j} '.mat'], 'file'))
                Quality = myload([direct img_names{j} '.mat'], 'Quality');
                if(~iscell(Quality)) % for backward data compatibility
                    Quality = {Quality};
                end
                
                thisQ{j} = Quality;
            else
                disp('didnt find a file');
                [direct img_names{j} '.mat']
            end
        end
        
        if(isempty(Q))
            Q = thisQ;
        else
            for j=1:numel(thisQ)
                for k = 1:numel(thisQ{j})
                    for l=1:numel(thisQ{j}{k})
                        Q{j}{k}(l).q = [Q{j}{k}(l).q; thisQ{j}{k}(l).q];
                        Q{j}{k}(l).sz = [Q{j}{k}(l).sz; thisQ{j}{k}(l).sz];
                    end
                end           
            end
        end
        %Q = cellfun(@conc, Q, thisQ, 'UniformOutput', false);
    end
    
    counter = 1;

    if(~isempty(class_label))
        cells_to_remove = false(numel(Q),1);
        for j=1:numel(Q) % iterate over images
            to_remove = ([Q{j}{1}(:).class] ~= class_label);
            
            Q{j}{1}(to_remove) = [];
            if(isempty(Q{j}{1}))
                cells_to_remove(j) = true;
            end
        end
        
        Q(cells_to_remove) = [];
    end
    
    avg_img_overlap = zeros(numel(Q),1);    
    for i=1:numel(Q) % iterate over images   
        avg_gt_overlap = zeros(numel(Q{i}),1);
        for j=1:numel(Q{i})  % iterate over gt segmentations
            max_gt_overlaps = zeros(numel(Q{i}{j}),1);
            mean_gt_overlaps = max_gt_overlaps;
            for k=1:numel(Q{i}{j}) % iterate over regions in the gt segmentation
                if(Q{i}{j}(k).class<255)
                    if(~isempty(Q{i}{j}(k).q))
                        m(counter) = max(Q{i}{j}(k).q);                    
                        avg(counter) = mean(Q{i}{j}(k).q);
                        mean_gt_overlaps(k) = mean(Q{i}{j}(k).q);
                    else
                        m(counter) = 0;
                        avg(counter) = 0;
                        mean_gt_overlaps(k) = 0;
                    end
                    
                    max_gt_overlaps(k) = m(counter);
                    counter = counter + 1;
                end
            end
            
            avg_gt_overlap(j) = mean(max_gt_overlaps);
        end
        avg_img_overlap(i) = mean(avg_gt_overlap);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['Dataset Avg Best ' segm_quality_type ':']);
    avg_best = mean(m)
    std_dev = std(m)
    avg_avg = mean(avg)
    
    avg_img_overlap(isnan(avg_img_overlap)) = [];

    disp(['Img Avg Best' segm_quality_type ':']);
    avg_best = mean(avg_img_overlap)
    std_dev = std(avg_img_overlap)
    
    disp(['Img Best Segmentation Covering with ' segm_quality_type ' score (avg over multiple GT):']);
    best_single_coverage = zeros(numel(Q), 1);
    best_single_coverage_best_GT = zeros(numel(Q), 1);
    for i=1:numel(Q) % iterate over images
        this_image_best_single_coverage = zeros(numel(Q{i}),1);
        for j=1:numel(Q{i}) % iterate over gt segmentations
            max_val = [];
            n_pixels_in_segments = [];
            total_n_pixels_gt = 0;
            the_k= [];
            for k=1:numel(Q{i}{j}) % iterate over regions in the image
                if(Q{i}{j}(k).class >= 255)
                    continue;
                end

                if(~isempty(Q{i}{j}(k).q))
                    [max_val(k), id] = max(Q{i}{j}(k).q);
                    n_pixels_in_segments(k) = Q{i}{j}(k).sz(id);
                else
                    max_val(k) = 0;
                    n_pixels_in_segments(k) = 0;
                end
                total_n_pixels_gt= total_n_pixels_gt + Q{i}{j}(k).object_sz;
                the_k = [the_k k];
            end
            max_val = max_val';
            
            if(~isempty(the_k))
                this_image_best_single_coverage(j) = (1/total_n_pixels_gt) * ([Q{i}{j}(the_k).object_sz]*max_val);
            else
                this_image_best_single_coverage(j) = -1;
            end
        end
        best_single_coverage(i) = mean(this_image_best_single_coverage); % average over humans
        
        if(isempty(this_image_best_single_coverage))
            this_image_best_single_coverage = -1;
        end
        best_single_coverage_best_GT(i) = max(this_image_best_single_coverage); % best over humans
    end
    
    % remove the ones we haven't computed
    best_single_coverage_best_GT(best_single_coverage_best_GT == -1) = [];
    best_single_coverage(best_single_coverage==-1) = [];
    best_single_coverage(isnan(best_single_coverage)) = [];
    avg_best_covering = mean(best_single_coverage)
    std_dev_covering = std(best_single_coverage)
    
    disp(['Img Best Segmentation Covering with ' segm_quality_type ' score (best over multiple GT) :']);
    avg_best_single_covering = mean(best_single_coverage_best_GT)
    std_dev_single_covering = mean(best_single_coverage_best_GT)
    
    n_segms = 0;
    num_Q = 0;
    for i=1:numel(Q) % count the number of segments we've got
        if(~isempty(Q{i}))
            n_segms = n_segms + numel(Q{i}{1}(1).sz);
            num_Q = num_Q + 1;
        end
    end
    
    avg_num_segments = sum(n_segms) / num_Q
    
    %visualize_bad_examples(exp_dir, segm_name{1}, img_names, best_single_coverage)

    
    %save([statistics_dir cell2mat(mask_type) '.mat'], 'avg_best', 'std_dev', 'avg_best_covering', 'std_dev_covering', ...
    %    'avg_num_segments');
end

function visualize_bad_examples(exp_dir, segm_name, img_names, best_single_coverage_best_GT, overlap_type)
    [sorted_vals, ids] = sort(best_single_coverage_best_GT, 'descend');
    %[sorted_vals, ids] = sort(best_single_coverage_best_GT);
        
    for i=1:numel(ids)
        I = imread([exp_dir 'JPEGImages/' img_names{ids(i)} '.jpg']);
        
        Quality = myload([exp_dir 'SegmentEval/' segm_name '/' overlap_type '/' img_names{ids(i)} '.mat'], 'Quality');
        if(~iscell(Quality)) % for backward data compatibility
            Quality = {Quality};
        end
        
        % load segments
        sorted_vals(i)
        load([exp_dir 'MySegmentsMat/' segm_name '/' img_names{ids(i)} '.mat']);        
        SvmSegm_show_best_segments(I,Quality, masks); figure;
        %subplot_auto_transparent(masks, I);
        
        % load all ground truth segmentations
        files = dir([exp_dir 'SegmentationObject/' img_names{ids(i)} '*']);
        %files = dir([exp_dir 'SegmentationSurfaces/' img_names{ids(i)} '*']);
        for j=1:numel(files)
            %filenames{j} = [exp_dir 'SegmentationSurfaces/' files(j).name];
            filenames{j} = [exp_dir 'SegmentationObject/' files(j).name];
            %gt{j} = imread([exp_dir 'SegmentationObject/' files(i).name]);
        end
        
        map = VOClabelcolormap(256);
        subplot(1,2,1), montage(filenames, map);
        subplot(1,2,2), sc(I);
        pause;
        close all;
    end
end

function res = conc(old,new)     
   res = old;
   if(isfield(old, 'q'))      
       for i=1:numel(new)
            res(i).q = [old(i).q; new(i).q];
            res(i).sz = [old(i).sz; new(i).sz];
       end
   else
       res = new;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% histogram gt segments with respect to their size, then see avg quality per bin
%  disp(['Segment Quality with gt segment size']);
%  duh = [Q{:}];
%  gt_sizes = [duh(:).object_sz];
%  [pos, sz] = hist(log10(gt_sizes), 30)
%  [pos, bin] = histc(log10(gt_sizes), [0 sz inf]);
%
%  10.^sz
%  for i=1:30
%    qual_hist(i) = sum(m(bin==i))/sum(bin==i);
%  end
%  qual_hist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
