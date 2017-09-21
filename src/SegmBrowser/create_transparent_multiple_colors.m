%function I = create_transparent_multiple_colors(segments, I, titles)
function merged_Img = create_transparent_multiple_colors(segments, I, use_voc_colors, labels, transparency, intensities) 
  DefaultVal('*use_voc_colors', 'false');
  DefaultVal('*transparency', '0.8');
  DefaultVal('*intensities', '[]');
  
  if(isempty(intensities))
      intensities = 0.5*ones(size(segments,3),1);
  end
  
%   if(size(segments,3) == 1 && numel(unique(segments)) > 1)
%       segments = labels2binary(segments);
%       segments(:,1) = []; % remove 0 label
%   end
  
  if(length(size(I))==2)
    I = repmat(I, [1 1 3]);
  end
  if(issparse(segments))
    segments = full(segments);
  end
  
  if(size(I,1) == size(segments,1))      
      segments = double(reshape(segments, [size(segments,1)*size(segments,2)], size(segments,3)));
  end
  
  n_segms = size(segments,2);
  
  merged_Img = I;
  bground = zeros(size(merged_Img,1), size(merged_Img,2), size(merged_Img,3));
 
  assert(size(segments,2) == numel(intensities));
    
  for i=1:n_segms
      if(use_voc_colors)
          if(exist('labels', 'var') && ~isempty(labels))
              if(labels(i) == 21) % background
                  bground = get_new_bg_voc(bground, segments(:,i), 0);
              else
                  bground = get_new_bg_voc(bground, segments(:,i), labels(i));
              end
          else
              bground = get_new_bg_voc(bground, segments(:,i), i);
          end
      else
        bground = get_new_bg_intensities(bground, segments(:,i), intensities(i));
      end
    
    alpha_chann = double(reshape(segments(:,i), size(I,1), size(I,2)))*transparency;
    merged_Img = immerge(merged_Img, bground, alpha_chann);
  end  
end

function bground = get_new_bg_voc(previous_bg, mask, i)
    bground = zeros(size(previous_bg));
    
    cmap = VOClabelcolormap(256);
    
    bground(:,:,1) = cmap(i+1,1)*255;
    bground(:,:,2) = cmap(i+1,2)*255;
    bground(:,:,3) = cmap(i+1,3)*255;
    
    mask_rs = repmat(reshape(mask,size(bground,1), size(bground,2)),[1 1 3]);
    bground = bground.*mask_rs;
    
    bground = previous_bg + bground;
end

function bground = get_new_bg_intensities(previous_bg, mask, intensity)    
    cmap = colormap('jet');
    intensity = max(1, round((intensity*size(cmap,1))));

    intensity = max(0, intensity);
    
    %end
    %max_intens = max(intensities);
    bground = zeros(size(previous_bg));
    bground(:,:,1) = cmap(intensity,1);
    bground(:,:,2) = cmap(intensity,2);
    bground(:,:,3) = cmap(intensity,3);
    bground = 255*bground;
    
    mask_rs = repmat(reshape(mask,size(bground,1), size(bground,2)),[1 1 3]);
    bground = bground.*mask_rs;
    
    bground = previous_bg + bground;
end

function bground = get_new_bg(previous_bg, mask, i)
  bground = zeros(size(previous_bg));
  if(i==1)
    bground(:,:,2) = 255;
  elseif(i==2)    
    bground(:,:,1) = 255;
  elseif(i==3)
    bground(:,:,3) = 255;
  elseif(i==4)
    bground(:,:,1) = 255;
    bground(:,:,2) = 255;
  elseif(i==5)
    bground(:,:,1) = 255;
    bground(:,:,3) = 255;
  elseif(i==6)
    bground(:,:,2) = 255;
    bground(:,:,3) = 255;    
  elseif(i==7)
    bground(:,:,1) = 125;
    bground(:,:,2) = 255;        
  elseif(i==7)
    bground(:,:,1) = 125;
    bground(:,:,3) = 255;            
  elseif(i==8)
    bground(:,:,2) = 125;
    bground(:,:,3) = 255;            
  elseif(i>=9)
    bground(:,:,1) = rand()*255;
    bground(:,:,2) = rand()*255;
    bground(:,:,3) = rand()*255;
  end
  
  mask_rs = repmat(reshape(mask,size(bground,1), size(bground,2)),[1 1 3]);
  bground = bground.*mask_rs;

  bground = previous_bg + bground;    
end