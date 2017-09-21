function error = hausdorff_surface_error_mex(vertices1, tri1, vertices2, tri2)    
    % there is normalization by the bounding volume diagonal of second model     
    r = 1000000000*rand(2,1) + numel(vertices1);
    str{1} = int2str(r(1));
    str{2} = int2str(r(2));
    
    assert(~strcmp(str{1}, str{2}));
    
    file1 = [str{1} '.off'];
    file2 = [str{2} '.off'];
    write_off(file1, vertices1, tri1);
    write_off(file2, vertices2, tri2);            
    
    error = meshx(file1,file2);
    
    delete(file1);
    delete(file2);
end
    
