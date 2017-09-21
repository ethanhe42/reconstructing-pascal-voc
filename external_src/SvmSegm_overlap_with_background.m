function SvmSegm_overlap_with_background(img_names, exp_dir, mask_type, segm_quality_type) %, opts)        
    segm_eval_dir = [exp_dir 'SegmentEval/'];

    if(~exist(segm_eval_dir))
        mkdir(segm_eval_dir);
    end

    direct = [segm_eval_dir mask_type '/' segm_quality_type '_withBG/' ];
    if ~exist(direct, 'dir')
        mkdir(direct);
    end

    t1 =  tic();
    %parfor i=1:length(img_names)
    for i=1:length(img_names)
        if(exist([direct img_names{i} '.mat'], 'file'))
            continue;
        end

        masks = myload([exp_dir 'MySegmentsMat/' mask_type '/' img_names{i} '.mat'], 'masks');

        name = [exp_dir 'SegmentationObject/' img_names{i} '.png'];
        alternatives = dir([name '*']);

        name =  [exp_dir 'SegmentationClass/' img_names{i} '.png'];
        class_alternatives = dir([name '*' ]);

        assert(~isempty(alternatives));
        assert(numel(alternatives)==numel(class_alternatives));
        

        ground_truth_obj_segs = cell(numel(alternatives),1);
        
        Quality = cell(1, numel(alternatives));
        for j=1:numel(alternatives)
        %parfor j=1:numel(alternatives)
            %Quality = struct();
            
            [classI, map] = imread([exp_dir 'SegmentationClass/' class_alternatives(j).name]);
            % repeated code, should substitute by function

            ground_truth_obj_segs{j} = imread([exp_dir 'SegmentationObject/'  alternatives(j).name]);

            un = unique(ground_truth_obj_segs{j});
            un(un~=0) = [];
            care = (classI~=255);
            for k=1:numel(un)
                ground_truth_k = zeros(size(ground_truth_obj_segs{j}));
                ground_truth_k(ground_truth_obj_segs{j} == un(k)) = 1;

                %disp('debugginnnnnn')
                %sz_masks = size(masks)
                %sz_gt = size(ground_truth_k)

                [duh1, duh2, qual] = myCalcCandScoreFigureGroundAll(masks,ground_truth_k, segm_quality_type, care);

                %%%% debugging %%%%
                %for m=1:200
                %    duh(m) = sum(masks(:,:,m) & ground_truth_k) / sum((masks(:,:,m) | ground_truth_k)) - qual(m);
                %end
                %mean(abs(duh))
                %max(abs(duh))
                
                new_un = unique(ground_truth_k);
                new_un(new_un==0) = [];
                clslbl = unique(classI(ground_truth_k==new_un)); % can do faster than this
                %assert(numel(clslbl) == 1);

                %img_names{i}
                %clslbl
                %return;
                Quality{j}(k).class = 21;
                Quality{j}(k).object = k;
                Quality{j}(k).object_sz = sum(sum(ground_truth_k));
                Quality{j}(k).image_name = img_names{i};

                if(isfield(Quality{j}(k), 'q') && ~isempty(Quality{j}(k).q))
                    Quality{j}(k).q
                    error('test this first');
                    for l=1:numel(Quality{j}(k).q)
                        if(Quality{j}(k).q(l) < qual(l))
                            Quality{j}(k).q(l) = qual(l);
                            Quality{j}(k).sz(l) = sum(sum(masks(:,:,l)));
                        end
                    end
                    Quality{j}(k).sz = sum(masks,3);
                else
                    Quality{j}(k).q = qual;
                    Quality{j}(k).sz = squeeze(sum(sum(masks,1), 2));
                end
            end
        end

        file_to_save = [direct img_names{i} '.mat'];
        mysave(file_to_save, 'Quality', Quality);
    end
    toc(t1)

end
