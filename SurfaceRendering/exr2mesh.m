height = exrread('output/regularized.exr');
height =  height(:,:,1);

y = 1:671;
x = 1:457;
[X,Y] = meshgrid(x,y);
s = surf(X,Y,height);
set(s,'LineStyle','none')
colorbar;
xlabel('X');
ylabel('Y');
zlabel('Z');

figure(2)
mesh(X,Y,height)
saveobjmesh('regularized_mesh.obj',X,Y,height)
