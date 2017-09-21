function the_var = myload(file, var)
    the_var = load(file, var);
    the_var = the_var.(var); % worst naming ever
end