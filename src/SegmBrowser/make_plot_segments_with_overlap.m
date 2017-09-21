function newI = make_plot_segments_with_overlap(masks, I, GT_img, voc_labels, the_cmap, line_width)
    DefaultVal('*the_cmap', '[]');
    DefaultVal('*GT_img', '[]');
    DefaultVal('*voc_labels', '[]');
    DefaultVal('*line_width', '2');
    
    if(isempty(the_cmap))
        cmap = VOClabelcolormap(256);
    else
        cmap = the_cmap;
        if(~all(the_cmap(1,:) ==[0 0 0]))
            the_cmap = [0 0 0; the_cmap];
        end
    end
    
    if(isempty(GT_img))
        GT_img = zeros(size(I,1), size(I,2));
    end
    
    if(iscell(voc_labels))
        VOCinit();
        for i=1:numel(voc_labels)
            new_voc_labels(i) = find(strcmp(VOCopts.classes, voc_labels{i}));
        end
        voc_labels = new_voc_labels;
    end
    
    if(~isempty(voc_labels))
        voc_labels(voc_labels==21) = 0;
        for i=1:numel(size(masks,3))
            GT_img(masks(:,:,i)) = voc_labels(i)+1;
        end
    end
    
    % get same color as for GT_IMG
    edge_I = zeros(size(masks));
    one_img_edge_I = zeros(size(GT_img));
    for i=1:size(masks,3)
       b{i} = cell2mat(bwboundaries(masks(:,:,i)));
       this_edge_I = zeros(size(GT_img));
       
       ind = sub2ind(size(GT_img), b{i}(:,1), b{i}(:,2));
       this_edge_I(ind) = 1;
       
       if(line_width==0)
           this_edge_I = logical(this_edge_I);
       elseif(line_width==1) % beautiful
           this_edge_I = bwmorph(this_edge_I, 'dilate');
       else
           this_edge_I = bwmorph(this_edge_I, 'dilate');
           this_edge_I = bwmorph(this_edge_I, 'dilate');
       end
       
       edge_I(:,:,i) = this_edge_I;
       
       one_img_edge_I(this_edge_I) = i+1;
    end
    %imshow(one_img_edge_I, cmap);
    
    newI = I;
       
    if(~isempty(voc_labels))
        one_img_edge_I(one_img_edge_I~=0) = 2;
    end        
    one_img_edge_I= ind2rgb(one_img_edge_I, cmap);
    one_img_edge_I = uint8(one_img_edge_I*255);   
    
    for i=1:3
        newI_i = newI(:,:,i);
        one_img_i = one_img_edge_I(:,:,i);
        to_change = sum(one_img_edge_I,3)~=0;
        newI_i(to_change) = one_img_i(to_change);       
        newI(:,:,i) = newI_i; 
    end
end