function [R,T,best_en] = cam_refinement(Rot,Tr,Shape,mask,projections,use_dt)
    % Refines the camera position based on the mask, so that all points
    % reproject inside the mask
    % minimizes:
    %weigths.proj*(keypoints - (rot*shape+tr))^2 + weights.dt*dist_mask(rot*shape +tr);

    if(use_dt)
        weights.proj = 1;
        weights.dt = 1;
    else
        weights.proj = 1;
        weights.dt = 0;
    end

    r = Rot(1:2,:);

    t = Tr(1:2);

    dt = get_distance_transform(mask,weights);

    %Precompute part of the gradient

    n_points = size(Shape,2);

    grad_t1 = [ones(n_points,1); zeros(n_points,1)];
    grad_t2 = [zeros(n_points,1); ones(n_points,1)];

    grad_r1 = [Shape(1,:) zeros(1,n_points)]';
    grad_r2 = [Shape(2,:) zeros(1,n_points)]';
    grad_r3 = [Shape(3,:) zeros(1,n_points)]';

    grad_r4 = [zeros(1,n_points) Shape(1,:)]';
    grad_r5 = [zeros(1,n_points) Shape(2,:)]';
    grad_r6 = [zeros(1,n_points) Shape(3,:)]';

    projections = projections(:);
    nan_proj = isnan(projections);


    step = 10^-4;
    n_iter=20000;

    coord = bsxfun(@plus,r*Shape,t);
    coord = coord';

    best_en=get_energy(coord,projections,dt,weights);
    en = best_en+1;

    dt_fx = zeros(n_points, 1);
    dt_fy = dt_fx;

    for k =1:n_iter
        coord = bsxfun(@plus,r*Shape,t);
        coord = coord';

        if(weights.dt~=0)
            [dt_fx,dt_fy]  = get_dt_gradient(dt,coord);
        end

        repr = weights.proj*2*(coord(:) - projections);
        repr(nan_proj) = 0;

        grad = [dt_fx;dt_fy] + repr;

        grad_r = [sum(grad.*grad_r1) sum(grad.*grad_r2) sum(grad.*grad_r3);...
            sum(grad.*grad_r4) sum(grad.*grad_r5) sum(grad.*grad_r6)];

        grad_t = [sum(grad.*grad_t1); sum(grad.*grad_t2)];

        while(en>best_en)
            if(step < 10^-7)
                if(any(grad_r(:)~=0))
                    step = 10^-2;
                    grad_r(:) = 0;
                else
                    break
                end
            end

            r1 = r-step*grad_r;
            r1 = projStiefel(r1);

            t1 = t-step*grad_t;
            c = bsxfun(@plus,r1*Shape,t1)';
            en=get_energy(c,projections,dt,weights);
            if(isnan(en)); en = best_en+1; end
            step = step/2;
        end

        if(en>best_en && step < 10^-7 && all(grad_r(:)==0))
            break;
        end
        
        step = step*4;
        r = r1;
        t = t1;

        best_en = en;

        if(best_en==0)
            break;
        end
        
        en = best_en+1;

    end

    scale = norm(r(1,:));
    R_3 = cross((1/scale)*r(1,:), (1/scale)*r(2,:));
    R = [r; scale*R_3];

    T =[t;0];


    end

    function en=get_energy(coord,projections,dt,weights)

    en = (coord(:) - projections).^2;
    en = weights.proj*sum(en(~isnan(en)));

    if(weights.dt ~=0)
        coord(:,2) = round(coord(:,2))+floor(dt.tolerance/2)*dt.M;
        coord(:,1) = round(coord(:,1))+floor(dt.tolerance/2)*dt.N;

        if(any(coord(:,1)<1) || any(coord(:,1)>size(dt.function,2)) || any(coord(:,2)>size(dt.function,1)) || any(coord(:,2)<1))
            en = nan;
            return;
        end

        a = dt.function(sub2ind(size(dt.function),coord(:,2),coord(:,1)));

        en = en + sum(a);
    end

    end


    function r = projStiefel(r)

    [U,D,V] = svd(r,'econ');
    c = (D(1,1) + D(2,2))/2;
    %c = mean(diag(D));

    r = c*U*V';

    end

    function [dt_fx,dt_fy]  = get_dt_gradient(dt,coord)


    ind = sub2ind(size(dt.function),round(coord(:,2))+floor(dt.tolerance/2)*dt.M,round(coord(:,1))+floor(dt.tolerance/2)*dt.N);

    dt_fx = dt.fx(ind);
    dt_fy = dt.fy(ind);
end

function distance_transform = get_distance_transform(mask,weights)
    tol = 3;

    distance_transform.tolerance = tol;

    [M,N] = size(mask);

    distance_transform.M = M;
    distance_transform.N = N;


    aux = zeros(tol*M,tol*N);
    aux(floor(tol/2)*M+1:ceil(tol/2)*M,floor(tol/2)*N+1:ceil(tol/2)*N) = mask;
    distance_transform.function = weights.dt*double(bwdist(aux)).^2;
    [distance_transform.x,distance_transform.y] = meshgrid(1:tol*N,1:tol*M);
    distance_transform.x = distance_transform.x - floor(tol/2)*N;
    distance_transform.y = distance_transform.y - floor(tol/2)*M;

    [distance_transform.fx,distance_transform.fy] = gradient(distance_transform.function);
end
