%function h = subplot_auto_transparent(segments, I, titles)
function [Imgs] = subplot_auto_transparent_parts(whole_segment, part_segments, I, titles, grid_type)
  if(isempty(part_segments))
      Imgs = [];
      return;
  end
  
  part_segments(part_segments==inf) = 10000;
  
  border_side = 0.05;
  border_top = 0.05;
  if(length(size(I))==2)
    I = repmat(I, [1 1 3]);
  end
  if(issparse(part_segments))
    part_segments = full(part_segments);
  end
  
  if(size(part_segments,1) ~= size(I,1))
    n_imgs = size(part_segments,2);
  else
    n_imgs = size(part_segments,3);
    part_segments = reshape(part_segments, size(part_segments,1) * size(part_segments,2), n_imgs);
  end

  n_rows = round(sqrt(n_imgs));
  n_cols = ceil(n_imgs/n_rows);

  assert(n_cols*n_rows >= n_imgs);

  counter = 1;
  
  % sets the color of the whole
  bground = zeros(size(I,1),size(I,2),3);
  bground(:,:,3) = 255;

  % sets the color of the parts
  bground_parts = zeros(size(I,1),size(I,2),3);
  bground_parts(:,:,1) = 255;
  bground_parts(:,:,2) = 255;

  
  Imgs = cell(n_imgs,1);
  
  
  if(size(part_segments,1) ~= size(I,1))
      alpha_chann_whole = reshape(whole_segment, size(I,1), size(I,2))*0.5;
  else
      alpha_chann_whole = whole_segment*0.5;
  end

  for i=1:n_imgs        
    if(size(part_segments,1) ~= size(I,1))
        alpha_chann_part = reshape(part_segments(:,i), size(I,1), size(I,2))*0.5;
    else
        alpha_chann_part = part_segments(:,:,i)*0.5;
    end
    
    Imgs{i} = immerge(I, bground, alpha_chann_whole);
    Imgs{i} = immerge(Imgs{i}, bground_parts, alpha_chann_part);

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