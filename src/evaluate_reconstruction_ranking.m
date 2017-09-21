% sometimes this crashes, (something related to parfor and writing and reading to files)
% just remove any .off files and rerun and it should be fine.
function evaluate_reconstruction_ranking(exp_dir,reconstr_name,withICP)                      
    if(nargin<3)
        withICP = 0;
    end
    
    for i=1:20
        classes{i} = VOC09_id_to_classname(i);
    end
    
    all_errors = cell(numel(classes),1);
    all_errors_upperbound = [];
    for i=1:20
              
       filename = ['./Results/ranking_errors_' classes{i} '_' reconstr_name '.mat'];
        
       if exist(filename, 'file')
           continue;
       end
       
       i       
       in_class = [];
        
       base_dir = ['./Results/' reconstr_name '/' classes{i} '/'];        
       load([base_dir 'reconstruction_data.mat']);
       
       sb = SegmBrowser(exp_dir, 'ground_truth', imgset_all);
       files = dir([base_dir '*.mat']);
       files(arrayfun(@(a) strcmp(a.name, 'reconstruction_data.mat'), files)) = [];

       t_class = tic();
       scores = [];
       class_errors = cell(numel(files),1);
       class_errors_ICP = class_errors;
       scores = class_errors;
       for j=1:numel(files)
           j
           gt_mesh = [];

           % load reconstructed mesh
           reconstructions = myload([base_dir files(j).name], 'reconstructions');
           img_name = sb.img_names{sb.whole_2_img_ids(in_class(reconstructions(1).triples(1)))};

           scores{j} = [reconstructions(:).score];

           % load ground truth synthetic mesh
           if 1
               load([exp_dir 'MyMeshes/ground_truth/' img_name '.mat'], 'tri', 'vertices');
           else
               tri = myload([exp_dir 'MyMeshes/ground_truth/' img_name '.mat'], 'tri');  % was trying to use parfor but couldn't
               vertices = myload([exp_dir 'MyMeshes/ground_truth/' img_name '.mat'], 'vertices');
           end

           gt_mesh.faces = tri';
           gt_mesh.vertices = vertices;

           errors = zeros(numel(reconstructions),1);
           errors_ICP = errors;
           parfor k=1:numel(reconstructions)
               if(numel(reconstructions(k).faces)<2)
                   errors(k) = inf;
               else
                   if(withICP)
                       [icpR,icpT] = icp(gt_mesh.vertices,reconstructions(k).vertices,'Matching','kDtree');
                       icpCorrected = bsxfun(@plus,icpR*reconstructions(k).vertices,icpT);
                       errors_ICP(k) = hausdorff_surface_error_mex(icpCorrected, reconstructions(k).faces, gt_mesh.vertices, gt_mesh.faces);
                   end
                   
                   errors(k) = hausdorff_surface_error_mex(reconstructions(k).vertices, reconstructions(k).faces, gt_mesh.vertices, gt_mesh.faces);
               end
           end
           class_errors{j} = errors;
           class_errors_ICP{j} = errors_ICP;
       end
       
       time_per_class = toc(t_class)
       
       errors = class_errors;
       if(withICP)
           errors_ICP = class_errors_ICP;
           save(filename, 'errors','errors_ICP','scores');
       else
           save(filename, 'errors','scores');
       end
    end   

    RANDOM_SCORE = [false true];
    %plotStyle = {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b--', 'g--', 'r--', 'c--', 'm--', 'y--', 'k--', 'b:', 'g:', 'r:', 'c:', 'm:', 'y:'};
    plotStyle = {'b', 'g', 'r', 'c', 'm', 'b', 'g', 'r', 'c', 'm', 'b', 'g', 'r', 'c', 'm', 'b', 'g', 'r', 'c', 'm'};
    
    all_errors_first = cell(2,20);    
    n_samples = [];
    for i=1:20
        subplot(4,5,i);        
        for k=1:numel(RANDOM_SCORE)            
            filename = ['./Results/ranking_errors_' classes{i} '_' reconstr_name '.mat'];        
            var = load(filename);  

            scores = var.scores;
            errors = var.errors;
            %errors = var.errors_ICP;
                
            n_samples = [n_samples; cellfun(@numel, scores)];
            
            % fill in any missing values
            n_proposals = max(cellfun(@numel, scores));
                        
            for l=1:numel(scores)
                if(numel(scores{l}) ~= n_proposals)
                    scores{l} = [scores{l} -inf(1, n_proposals - numel(scores{l}))];
                    errors{l} = [errors{l}; max(errors{l})*ones(n_proposals - numel(errors{l}),1)];
                end
            end
            
            scores = cell2mat(scores)';
            errors = cell2mat(errors');

            errors(errors==inf) = 100;

            if RANDOM_SCORE(k)
                scores = rand(size(scores));
            end

            if(size(errors,1)>1)
                [~, ids_srt] = sort(scores, 'descend');
                for j=1:size(errors,2)
                    errors(:,j) = errors(ids_srt(:,j),j);
                end
            end
            
            % get min of avg error per rank
            
            min_avg = [];
            %if(fix_me)
            %    errors = errors(1,:);
            %end
            
            for j=1:size(errors,1)
                if(size(errors,1)>1)
                    min_avg(j) = mean(min(errors(1:j,:),[],1));
                else
                    min_avg(j) = mean(errors);
                end

                if j==1
                    all_errors_first{k,i} = min(errors(1:j,:),[],1);
                end
            end
            
            if(size(errors,1)>1)
                all_errors_upperbound{i} = min(errors);
            else
                all_errors_upperbound{i} = errors;
            end
            
            if RANDOM_SCORE(k)
                h(i) = plot(1:numel(min_avg), min_avg, [plotStyle{i} 'o'], 'LineWidth', 3); hold on;     
            else
                h(i) = plot(1:numel(min_avg), min_avg, plotStyle{i}, 'LineWidth', 5); hold on;  
            end                                    
        end
        legend({[classes{i} ' - predicted'], [classes{i} ' - random']});
    end
    
    figure;
    cmap = VOClabelcolormap(256);
    for i=1:20
        plot(i, mean(all_errors_first{2,i}), 'x', 'MarkerSize', 15, 'LineWidth', 5, 'Color', cmap(i+1,:)); hold on;
        plot(i, mean(all_errors_first{1,i}), 'o', 'MarkerSize', 15, 'LineWidth', 5, 'Color', cmap(i+1,:)); hold on;
        plot(i, mean(all_errors_upperbound{i}), 's', 'MarkerSize', 15, 'LineWidth', 5, 'Color', cmap(i+1,:)); hold on;
        fprintf('%s:\t\t %f \t %f \t %f\n', classes{i}, mean(all_errors_first{1,i}), mean(all_errors_upperbound{i}), mean(all_errors_first{2,i}));
        classes{i} = VOC09_id_to_classname(i);
    end    
    set(gca, 'Xlim', [0 21]);
    xticklabel_rotate(1:20, 45, classes, 'FontSize', 15)
    set(gca,'FontSize',15);
    legend('Random', 'Top-ranked', 'Best available');
    
    all_errors_upperbound = cell2mat(all_errors_upperbound);
    
    fprintf('\n\n');
    mean_error_pred = mean(cell2mat(all_errors_first(1,:))');
    mean_error_rand = mean(cell2mat(all_errors_first(2,:))');
    mean_error_upperbound = mean(all_errors_upperbound);
    fprintf('Mean error using predicted: %f\nMean error with random selection: %f\n', mean_error_pred, mean_error_rand);
    fprintf('Lower bound on error of a selection function: %f\n', mean_error_upperbound);
    fprintf('Average number of proposals per image: %f\n', mean(n_samples));
end