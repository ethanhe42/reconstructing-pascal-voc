function mysave(file, name, var)
    eval([name ' =  var;']);
    save(file, name, '-V6');   
end
