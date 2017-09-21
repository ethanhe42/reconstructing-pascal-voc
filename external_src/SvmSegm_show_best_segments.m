%function SvmSegm_show_best_segments(I, Q, masks)
% I is the image
% Q is the qualities of each segment (find the file in SegmentEval folder)
  % masks are the computed segments
function [best, max_score, best_seg_id] = SvmSegm_show_best_segments(I, Q, masks, n)
  DefaultVal('*n', '1');
  if(iscell(Q))
    Q = Q{1};
  end  
  for i=1:numel(Q)
    if(n~=1)
      [max_score{i}, best_seg_id{i}] = sort(Q(i).q, 'descend');
      best_seg_id{i} = best_seg_id{i}(1:n);
      tit = [];
    else    
      [max_score(i), best_seg_id(i)] = max(Q(i).q);
      tit{i} = sprintf('%f', Q(i).q(best_seg_id(i)));
    end
  end
  if(iscell(best_seg_id))
    best_seg_id = cell2mat(best_seg_id);
    max_score = cell2mat(max_score);
  end
  best = subplot_auto_transparent(masks(:,:,best_seg_id), I, tit);
  %figure;
  %sc(seg, 'rand')
end