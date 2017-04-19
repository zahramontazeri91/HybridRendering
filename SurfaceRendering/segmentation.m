EXRfile = 'input/height.exr';

height = exrread(EXRfile);
grayImage = height(:,:,1);
[x_size y_size] = size(grayImage);
imagesc(grayImage);
line([50 50], [1 670])
line([290 290], [1 670])
line([134 134], [1 670])
line([227 227], [1 670]) 
line([360 360], [1 670]) 
line([440 440], [1 670])

line([1 457], [115 115])
line([1 457], [220 220])
line([1 457], [335 335])
line([1 457], [440 440])
line([1 457], [565 565])

for i=12:12
    % Draw ROI
    [binaryMask, x, y] = roipoly();
    hold on;
    plot(x, y, 'r.-', 'MarkerSize', 15);
    
    filename = ['input/manually mask/patch_' num2str(i) '.png'];
    imwrite(binaryMask,filename);
end
