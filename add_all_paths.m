function exp_dir = add_all_paths()        
    
    addpath('./src/');
    addpath(['./src/SegmBrowser/']);

    addpath('./external_src/');
    addpath('./external_src/marques_costeira/');
    addpath('./external_src/RenderMe/Matlab/');
    addpath('./external_src/point_distance_triangle/');    
    addpath('./external_src/Mexsh-1.13/bin/');
    addpath('./external_src/VOCcode/');
    addpath('./external_src/immerge/');
    
    exp_dir = './VOC/';
    
    try
      rng(1234);
    catch
      disp('don''t have rng');
    end
    
    
    
    
end
