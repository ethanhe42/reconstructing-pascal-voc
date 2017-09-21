function [ap, fp, tp, ovmax] = voc_ap_fast_and_stupid(boxes1, class_bbox_overlaps)
  n_pred_per_img = cellfun(@(a) size(a,1), boxes1);
  cs = cumsum(n_pred_per_img);
  first = 1;
  
  n = numel(boxes1);
  ids = cell(n,1);
  for j=1:n
    if(size(boxes1{j},1)~=0)
      ids{j} = first:cs(j);
      first = cs(j)+1;
    else
      ids{j} = [];
    end
  end
  n_pred = sum(n_pred_per_img);
  fp = zeros(n_pred,1);
  tp = zeros(n_pred,1);

  all_of_them = cell2mat(class_bbox_overlaps');  
  npos = numel(all_of_them) - sum([all_of_them.diff]);
  
  % assign all predictions in images that do not have objects as
  % false positives
  non_obj_cells = cellfun(@(a) isempty(a), class_bbox_overlaps);
  fp(cell2mat(ids(non_obj_cells)')) = true;

  % assign each segment to the most overlapping object
  obj_cells = find(~non_obj_cells);
  for j=1:numel(obj_cells) % search only for images with gt objects of this class
    gt_taken = false(numel(class_bbox_overlaps{obj_cells(j)}),1);  
    ov_max = zeros(numel(ids{obj_cells(j)}),1);
    for k=1:numel(ids{obj_cells(j)}) % search over detections
      ov_max(k) = -inf;
      d = ids{obj_cells(j)}(k);
      
      for m=1:numel(class_bbox_overlaps{obj_cells(j)}) % seach over objects
        [this_overlap] = class_bbox_overlaps{obj_cells(j)}(m).q(k);
        if(this_overlap>ov_max(k))
          ov_max(k) = this_overlap;
          obj_max = m;
        end
      end

      if ov_max(k)>=0.5
        if ~class_bbox_overlaps{obj_cells(j)}(obj_max).diff
          if ~gt_taken(obj_max)
            tp(d)=1;% true positive
            gt_taken(obj_max)=true;
          else
            fp(d)=1;            % false positive (multiple detection)
          end
        end
      else
        fp(d)=1;                    % false positive
      end
    end
  end  

  bb_array = cell2mat(boxes1)';
  confidence = bb_array(5,:);
  [sc,si]=sort(-confidence);
    
  sorted_fp = fp(si);
  sorted_tp = tp(si);
  
  sum_fp=cumsum(sorted_fp);
  sum_tp=cumsum(sorted_tp);
  rec=sum_tp/npos;
  %prec=sorted_tp./(sum_fp+sum_tp);
  prec=sum_tp./(sum_fp+sum_tp);
  ap=VOCap(rec,prec);
end