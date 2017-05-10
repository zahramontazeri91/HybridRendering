clc
clear

fileID = fopen ('Output/Data/blockID_heightmap.dat', 'w');
blockID = [4,0,0,0,0];
fprintf(fileID, '%i', blockID);
fclose(fileID);

type Output/Data/blockID_heightmap.dat

%%
% read from file
fileID_vol = fopen ('Output/Data/gabardine_400x600_tiled2.dat');
%blockID_vol = fscanf(fileID_vol, '%f' );
blockID_vol = fscanf(fileID_vol, '%f' )
fclose(fileID_vol);