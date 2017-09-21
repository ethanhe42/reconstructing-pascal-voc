% returns bounding box even is mask is composed of multiple connected
% components
function rp_bbox = regionprops_BB_mine(mask, slack)
  DefaultVal('*slack', '0');
  
  [c_x, c_y] = find(mask);
  if(isempty(c_x)) % robust to masks with zero pixels.
    rp_bbox = [1 1 1 1];
  else
    rp_bbox(1) = min(c_y);
    rp_bbox(2) = min(c_x);
    rp_bbox(3) = max(c_y) - rp_bbox(1) + 1;
    rp_bbox(4) = max(c_x) - rp_bbox(2) + 1;
    
    if(slack~=0)      
      rp_bbox(1) = max(1, rp_bbox(1)-slack);
      rp_bbox(2) = max(1, rp_bbox(2)-slack);
      rp_bbox(3) = min(size(mask,2)-rp_bbox(1), rp_bbox(3)+2*slack);
      rp_bbox(4) = min(size(mask,1)-rp_bbox(2), rp_bbox(4)+2*slack);
    end
    
    % showboxes(mask,[rp_bbox(1) rp_bbox(2) rp_bbox(1)+rp_bbox(3) rp_bbox(2)+rp_bbox(4)])
  end
end