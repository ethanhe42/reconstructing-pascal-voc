function compile_all()
  
  % Compile the rendering mex file (requires OSMesa)
  cd external_src/RenderMe/Matlab/
  compile % this may produce some warnings
  cd ../../..
  
  % Compile the evaluation mex file 
  cd external_src/Mexsh-1.13/
  compile
  cd ../../
  
end

