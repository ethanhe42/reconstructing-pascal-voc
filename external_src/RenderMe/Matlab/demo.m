%{

This code is to render a Mesh given a 3x4 camera matrix with an image resolution widthxheight. The rendering result is an ID map for facets, edges and vertices. This can usually used for occlusion testing in texture mapping a model from an image, such as the texture mapping in the following two papers.

--Jianxiong Xiao http://mit.edu/jxiao/

Citation:

[1] J. Xiao, T. Fang, P. Zhao, M. Lhuillier, and L. Quan
Image-based Street-side City Modeling
ACM Transaction on Graphics (TOG), Volume 28, Number 5
Proceedings of ACM SIGGRAPH Asia 2009

[2] J. Xiao, T. Fang, P. Tan, P. Zhao, E. Ofek, and L. Quan
Image-based Facade Modeling
ACM Transaction on Graphics (TOG), Volume 27, Number 5
Proceedings of ACM SIGGRAPH Asia 2008

%}


clear
clc
close all

%compile
data()

t = DelaunayTri(vertex(1,:)', vertex(2,:)', vertex(3,:)');
face = uint32(t.Triangulation');
vertex = t.X';

%result = RenderMex2(P, img_width, img_height, vertex(:,1:1000), edge, face-1)';
%result = RenderMex(P, img_width, img_height, vertex, edge, face-1)';
result = RenderTriMex(P, img_width, img_height, vertex, edge, face(1:3,:)-1)';

close all
imagesc(result)
axis equal
axis tight
max(max(result))
