% for regularized map
regularized = exrread('Output/regularized.exr');
regularized =  regularized(:,:,1);
[m,n] = size(regularized);
y = 1:m;
x = 1:n;
%regularized = medfilt2(regularized, [7 7]);
regularized = imgaussfilt(regularized, 3);

[X,Y] = meshgrid(x,y);

% figure(1)
% axis equal
% mesh(X,Y,regularized)
% colorbar;
% xlabel('X');
% ylabel('Y');
% zlabel('Z');

saveobjmesh('output/regularized_mesh.obj',X,Y,regularized)

%%
% visualize residual map
% residual = exrread('Output/residual.exr');
% residual =  residual(:,:,1);
% 
% figure(2)
% imagesc(residual)
% colorbar;

%%

% test = imread('input/manually mask/patch_1.png');
% 	for i=2:30
%         filename = ['input/manually mask/patch_' num2str(i) '.png'];
% 		test = test + imread(filename);
% end
% imagesc(test);
