%function h = subplot_auto_transparent(segments, I, titles)
function [Imgs] = subplot_auto_transparent_parts_noov(whole_segment, part_segments, I, voc_color_ids, titles, grid_type)
  DefaultVal('*voc_color_ids', '[]');
  if(isempty(part_segments))
      Imgs = [];
      return;
  end
  
  if(~iscell(whole_segment))
      whole_segment = {whole_segment};
      part_segments  = {part_segments};
      I = {I};
  else
      assert(iscell(part_segments));
      assert(iscell(I));
  end
  
  border_side = 0.05;
  border_top = 0.05;

  n_rows = 1;
  n_cols = numel(I);
  counter = 1;  

  cmap = VOClabelcolormap(256);
  
  n_imgs = numel(I);
  Imgs = cell(n_imgs,1);
  for h=1:n_imgs
      if(length(size(I{h}))==2)
          I{h} = repmat(I{h}, [1 1 3]);
      end
      if(issparse(part_segments{h}))
          part_segments{h} = full(part_segments{h});
      end

      if(size(part_segments,1) ~= size(I,1))
          alpha_chann_whole = reshape(whole_segment{h}, size(I{h},1), size(I{h},2))*0.5;
      else
          alpha_chann_whole = whole_segment{h}*0.5;
      end
  
      Imgs{h} = I{h};
      
      bground = zeros(size(I{h},1),size(I{h},2),3);
      bground(:,:,3) = 255;
      
      if(isempty(part_segments{h}))
          Imgs{h} = immerge(Imgs{h}, bground, alpha_chann_whole);
      else
          for i=1:size(part_segments{h},3)
              alpha_chann_part = part_segments{h}(:,:,i)*0.7;
              
              if(~isempty(voc_color_ids))
                c_id = voc_color_ids{h}(i);
              else
                c_id = i;
              end
              bground_parts = zeros(size(I{h},1),size(I{h},2),3);
              bground_parts(:,:,1) = cmap(c_id, 1)*255;
              bground_parts(:,:,2) = cmap(c_id, 2)*255;
              bground_parts(:,:,3) = cmap(c_id, 3)*255;
              
              Imgs{h} = immerge(Imgs{h}, bground_parts, alpha_chann_part); %
              
              counter = counter + 1;
          end
      end
  end
  
  %montage_new(Imgs, titles, 'Size', [n_rows n_cols], 'Border', 0.1);
  if(~exist('grid_type', 'var'))
      grid_type = [];
  end
  
  if(exist('titles', 'var') && ~isempty(titles))
    if(~iscell(titles)) % if not a cell assumes they're numbers
        for i=1:length(titles)
            new_titles{i} = sprintf('%f', titles(i));
        end
        titles = new_titles;
    end
    if(~exist('grid_type', 'var'))
        montage_new(Imgs, titles, 'Border',  [border_top border_side]);
    else
        montage_new(Imgs, titles, 'Border', [border_top border_side], 'Size', grid_type);
    end
  else
    montage_new(Imgs, [], 'Border',  [border_top border_side], 'Size', grid_type);
  end
%   hold on;
%   for i=1:n_imgs
%     if(nargin==3)
%       h(counter) = subplot(n_rows,n_cols,i); title(titles{i});
%     end
%   end