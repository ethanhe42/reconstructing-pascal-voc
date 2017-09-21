function [newR,newT] = refine_all_cameras(Shape,R,T,masks,kp)

    n_images = size(R,3);

    newR = zeros(size(R));
    newT = zeros(size(T));

    for h=1:n_images
        thisR = R(:,:,h);
        thisT = T(:,h);
        mask = masks{h};
        thisKP = kp(:,:,h)';

        if(sum(~isnan(thisKP(:,1)))<=1)
            %There is one car image that only has one keypoint for which the
            %rotation estimation is just wrong
            continue;
        end


        %% refinement only based on the keypoints
        %It considers two different initializations, and chooses the solution with
        %best energy. The initialization coming from the rigid reconstruction and
        %identity rotation and translation given by the centroid of the mask.

        %1st initialization
        [refR,refT,en] = cam_refinement(thisR,thisT,Shape,mask,thisKP,0);


        centroid = regionprops(mask,'centroid');
        centroid = [centroid(1).Centroid';0];
        %2nd initialization
        [refR1,refT1,en1] = cam_refinement(eye(3),centroid,Shape,mask,thisKP,0);

        if(en<en1)
            thisR = refR;
            thisT = refT;
        else
            thisR = refR1;
            thisT = refT1;
        end

        %% Refinement with silhouette information
        [thisR,thisT] = cam_refinement(thisR,thisT,Shape,mask,thisKP,1);

        newR(:,:,h) = thisR;
        newT(:,h) = thisT;
    end
end
