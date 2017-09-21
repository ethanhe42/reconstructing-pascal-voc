%function h = subplot_auto_transparent(segments, I, titles)
function [Imgs] = subplot_auto_transparent_parts_multiple_imgs(whole_segments, part_segments, I, titles) 
  border_side = 0.05;
  border_top = 0.05;

  
  Imgs = I;
  for i=1:numel(I)
      if(length(size(I{i}))==2)
        I{i} = repmat(I{i}, [1 1 3]);
      end

      % sets the color of the whole
      bground = zeros(size(I{i},1),size(I{i},2),3);
      bground(:,:,3) = 255;
      
      % sets the color of the parts
      bground_parts = zeros(size(I{i},1),size(I{i},2),3);
      bground_parts(:,:,1) = 255;
      bground_parts(:,:,2) = 255;
      

      alpha_chann_whole = reshape(whole_segments{i}, size(I{i},1), size(I{i},2))*0.5;
      
      alpha_chann_part = part_segments{i}*0.5;
      
      Imgs{i} = immerge(I{i}, bground, alpha_chann_whole);
      Imgs{i} = immerge(Imgs{i}, bground_parts, alpha_chann_part);      
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
      end
  else
      montage_new(Imgs, [], 'Border',  [border_top border_side]);
  end
  %   hold on;
  %   for i=1:n_imgs
  %     if(nargin==3)
  %       h(counter) = subplot(n_rows,n_cols,i); title(titles{i});
  %     end
  %   end
end
  