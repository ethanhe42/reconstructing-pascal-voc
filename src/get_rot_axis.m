function [ ax ] = get_rot_axis(Shape,bot_top)

axes = {'x','y','z'};


a = mean(Shape(:,bot_top{1}),2) - mean(Shape(:,bot_top{2}),2);

[~,b] = max(abs(a));

ax = axes{b};


end

