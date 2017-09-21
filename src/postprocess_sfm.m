function [Shape, R, k] = postprocess_sfm(Shape, R)
    % scale Shape (set max value to k)
    m = max(abs(Shape(:)));
    k = 20;
    Shape = (k/m) * Shape;
    R =  (m/k) * R;
    
    if 1
        % Align Shape with axis                
        [P,duh, e] = princomp(Shape');
        
        P = P';
        if(abs(det(P)-1)>0.1)
            P(3,:) = -P(3,:);
        end
        
        Shape = P*Shape;
        for i=1:size(R,3)
            thisR = R(:,:,i);
            s = norm(thisR(1,:));
            thisR = thisR/s;
            thisR = thisR*P';
            if(any(isnan(thisR)))
                thisR = zeros(size(thisR));
            end
            R(:,:,i) = s*thisR;             
        end
    end    
end