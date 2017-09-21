function [scores] = get_ranking(not_projected, dist_proj)

un = unique(not_projected);

if(size(un,1)>1)
    not_projected(not_projected ==un(1)) = un(2);
    
    ranking = [];
    
    for k = 2:length(un);
        ind = find(not_projected == un(k));
        
        [a,b] = sort(dist_proj(ind));
        
        ranking = [ranking; ind(b)];
        
    end
    
else
    [~,ranking] = sort(dist_proj);
    
end

for i=1:numel(ranking)
    scores(i) = -find(ranking==i);
end

end