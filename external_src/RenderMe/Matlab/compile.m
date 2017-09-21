% one of these two should work
try
	mex -O -lOSMesa RenderTriMex.cpp 
catch
	mex -O -lOSMesa -lGL RenderTriMex.cpp 
end
