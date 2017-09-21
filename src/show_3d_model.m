function show_3d_model(Shape, kp_names, type, rigidity)    
    DefaultVal('*type', '[]');

    n_kp = size(Shape,2);
    colormap('hsv');
    %%%%%%%%%%%%% VISUALIZATION of 3d model %%%%%%%%%%%%%%
    if strcmp(type, 'convex_hull')
        %%%% visualize convex hull %%%
        try
          T = delaunay(Shape(1,:), Shape(2,:), Shape(3,:));
        catch
          T = delaunay3(Shape(1,:), Shape(2,:), Shape(3,:));
        end
        tetramesh(T, Shape'); %camorbit(20,0)
    elseif strcmp(type, 'rigid_links')
        %%%% visualize rigid links
        [a,b] = find(triu(rigidity,1));

        for i=1:numel(a)
            line([Shape(1,a(i)) Shape(1,b(i))], [Shape(2,a(i)) Shape(2,b(i))], [Shape(3,a(i)) Shape(3,b(i))]); hold on;
        end
        axis equal;
    else
        % visualize keypoint 3d locations in isolation
        plot3(Shape(1,:), Shape(2,:), Shape(3,:),'o', 'MarkerSize', 12, 'LineWidth', 8); hold on;
    end    

    if~isempty(kp_names)
        hold on;
        for i=1:n_kp
            h = text(Shape(1,i), Shape(2,i), Shape(3,i), kp_names{i}, 'FontSize', 10, 'Color', [0 0 0], 'BackgroundColor', [1 1 1]);
            set(h,'interpreter','none');
        end
    end
end