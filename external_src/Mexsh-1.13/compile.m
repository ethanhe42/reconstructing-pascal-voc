% make till it fails to compile everything, then mex to link it properly

system('make clean');
system(['make all MATLABDIR=' matlabroot]);

try
    mex -O obj/compute_error.o obj/model_analysis.o obj/mesh_run.o obj/reporting.o obj/colormap.o obj/xalloc.o obj/mesh.o obj/ScreenWidget.o obj/TextWidget.o obj/Error3DViewerWidget.o obj/ColorMapWidget.o obj/Lighted3DViewerWidget.o obj/Basic3DViewerWidget.o obj/InitWidget.o obj/moc_Basic3DViewerWidget.o obj/moc_Lighted3DViewerWidget.o obj/moc_Error3DViewerWidget.o obj/moc_ScreenWidget.o obj/moc_InitWidget.o obj/moc_ColorMapWidget.o lib/lib3d.a -L/usr/lib64/qt-3.3/lib -L/usr/X11R6/lib -lqt -lGL -lGLU -lpthread -lXmu -lXext -lX11 -lm -lz   -o bin/meshx
catch
    mex -O obj/compute_error.o obj/model_analysis.o obj/mesh_run.o obj/reporting.o obj/colormap.o obj/xalloc.o obj/mesh.o obj/ScreenWidget.o obj/TextWidget.o obj/Error3DViewerWidget.o obj/ColorMapWidget.o obj/Lighted3DViewerWidget.o obj/Basic3DViewerWidget.o obj/InitWidget.o obj/moc_Basic3DViewerWidget.o obj/moc_Lighted3DViewerWidget.o obj/moc_Error3DViewerWidget.o obj/moc_ScreenWidget.o obj/moc_InitWidget.o obj/moc_ColorMapWidget.o lib/lib3d.a -L/usr/lib64/qt-3.3/lib -L/usr/X11R6/lib -lqt-mt -lGL -lGLU -lpthread -lXmu -lXext -lX11 -lm -lz   -o bin/meshx    
end