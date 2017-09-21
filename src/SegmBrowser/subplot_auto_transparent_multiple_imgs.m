%function h = subplot_auto_transparent(segments, I, titles)
function [Imgs, handles] = subplot_auto_transparent_multiple_imgs(segments, I, titles, voc_cmap_ids, grid_type) 
  border_side = 0.05;
  border_top = 0.05;

  for i=1:numel(I)
      if(length(size(I{i}))==2)
        I{i} = repmat(I{i}, [1 1 3]);
      end

      counter = 1;

      if(exist('voc_cmap_ids', 'var') && ~isempty(voc_cmap_ids))      
          cmap = VOClabelcolormap();
          assert(size(segments,2) == numel(voc_cmap_ids));

          bground{i} = zeros(size(I,1),size(I,2),3);
          bground{i}(:,:,1) = cmap(voc_cmap_ids(i),1);
          bground{i}(:,:,2) = cmap(voc_cmap_ids(i),2);
          bground{i}(:,:,3) = cmap(voc_cmap_ids(i),3);
          bground{i} = uint8(255*bground{i});
      else                
          bground{i} = zeros(size(I{i},1),size(I{i},2),3);
          bground{i}(:,:,2) = 255;          
      end
  end
  
  Imgs = I;
  n_imgs = numel(Imgs);
  
  for i=1:n_imgs
      alpha_chann = segments{i}*0.5;
      
      %sc(sc(I).*sc(alpha_chann))    
      Imgs{i} = immerge(Imgs{i}, bground{i}, alpha_chann);
      
      counter = counter + 1;
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
        handles = montage_new(Imgs, titles, 'Border',  [border_top border_side]);
    else
        handles = montage_new(Imgs, titles, 'Border', [border_top border_side], 'Size', grid_type);
    end
  else
    handles = montage_new(Imgs, [], 'Border',  [border_top border_side], 'Size', grid_type);
  end
%   hold on;
%   for i=1:n_imgs
%     if(nargin==3)
%       h(counter) = subplot(n_rows,n_cols,i); title(titles{i});
%     end
%   end