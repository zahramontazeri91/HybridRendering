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
% visualize the final results
% regularized = exrread('Output/regularized.exr');
% regularized =  regularized(:,:,1);
% height = exrread('input/height.exr');
% height =  height(:,:,1);
% residual = exrread('Output/residual.exr');
% residual =  residual(:,:,1);
% 
% figure(2)
% subplot(1,3,1), imagesc(residual), title('Residual Map');
% caxis([0,170]), colorbar;
% subplot(1,3,2), imagesc(regularized), title('Regularized Map');
% caxis([0,170]), colorbar;
% subplot(1,3,3), imagesc(height), title('Height Map');
% caxis([0,170]), colorbar;
%%

% test = imread('input/auto manually mask/patch_1.png');
% 	for i=2:30
%         filename = ['input/auto manually mask/patch_' num2str(i) '.png'];
% 		test = test + imread(filename);
% end
% imagesc(test);

%%
    %kronecker tenser product for tiling:
    A = [1,1,1,1; 1,1,1,1];
    reg_height_tiled = kron(A, reg_height);
