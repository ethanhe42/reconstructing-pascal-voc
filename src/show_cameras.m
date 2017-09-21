function show_cameras(Shape, kp_names, R, T, I, masks)
    function myupdatefcn_simple(eventdata,duh)
      figure;
      id = get(eventdata, 'UserData');      
      subplot_auto_transparent(masks{id}, I{id});
    end

    if 1
        show_3d_model(Shape, [], 'convex_hull'); hold on;

        radius = 50;
        k = 1;

        x = k*[0 0 0 0; 1 1 -1 1; 1 -1 -1 -1];
        y = k*[0 0 0 0; 1 1 -1 -1; -1 1 1 -1];
        z = k*[0 0 0 0; 5 5 5 5; 5 5 5 5];

        point = [0 0 -radius]';
        z =  -radius + z;

        centers = zeros(3, size(R,3));    
        for h=1:size(R,3)
            thisR = R(:,:,h);

            s = norm(thisR(1,:));
            thisR = thisR/s;
            invR = thisR';
            
            
            this_x = zeros(size(x));
            this_y = this_x;
            this_z = this_x;
            for i=1:size(x,1)
                for j=1:size(x,2)
                    v = invR*[x(i,j) y(i,j) z(i,j)]';
                    %v = [x(i,j) y(i,j) z(i,j)]' - T(:,h);
                    %v = R(:,:,h)\v;
                    %v = bsxfun(@plus, R(:,:,h)\[x(i,j) y(i,j) z(i,j)]', T(:,h));
                    this_x(i,j) = v(1);
                    this_y(i,j) = v(2);
                    this_z(i,j) = v(3);
                end
            end
            centers(:,h) = [this_x(1) this_y(1) this_z(1)]';
            thish = fill3(this_x,this_y,this_z, rand(3,4)); hold on;
            set(thish, 'UserData', h, 'ButtonDownFcn',@myupdatefcn_simple);
            p1 = invR*point;
            p2 = invR*(-point);
            line([p1(1);p2(1)],[p1(2);p2(2)],[p1(3);p2(3)]);
        end

        %fill3(x,y,z, rand(3,4)); grid on;
        %view(100, 30)
        axis tight;
        axis equal;
        view(90,0);

        %xlim([-1000 1000]);
        %ylim([-1000 1000]);
        %zlim([-1000 1000]);
        grid on;

        %dcm_obj = datacursormode(fig);
        %set(dcm_obj,'UpdateFcn',@(a,b) myupdatefcn(a,b,sb,obj_ids, centers, dcm_obj));
    else
        for j=1:size(R,3)
            s(j) = norm(R(1,:,j));
            R(:,:,j) = R(:,:,j)/s(j);   
            R_inv(:,:,j) = inv(R(:,:,j));
        end
        
        anim = Animation();
        anim.K = [s; zeros(1,numel(s)); s];
        anim.S = Shape;
        anim.R = R_inv;
        anim.t = -T;
        playAnim(anim, 'frame', 1, 'nCam', 10 );
    end
end
