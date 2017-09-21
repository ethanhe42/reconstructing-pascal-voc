function init_rnd_seed()
    try % initialize random seed
        rng(1234);
    catch
        disp('don''t have rng');
    end
end